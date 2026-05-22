use std::process::{Command, Stdio};
fn ours() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_rsomics-peak-count"))
}
fn golden(n: &str) -> String {
    format!("{}/tests/golden/{}", env!("CARGO_MANIFEST_DIR"), n)
}

fn have(tool: &str) -> bool {
    Command::new(tool)
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok_and(|s| s.success())
}

#[test]
fn counts_are_non_negative() {
    let out = Command::new(ours())
        .arg(golden("small.bam"))
        .args(["-b", &golden("peaks.bed")])
        .output()
        .unwrap();
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    for line in s.lines() {
        let count: u64 = line.split('\t').nth(3).unwrap().parse().unwrap();
        assert!(count < 1_000_000, "count should be reasonable: {count}");
    }
}

// Per-peak read counts must match `bedtools multicov` (a read counts toward a
// peak it overlaps by >=1bp via its reference span).
#[test]
fn counts_match_bedtools_multicov() {
    if !have("bedtools") || !have("samtools") {
        eprintln!("skipping: bedtools/samtools not found");
        return;
    }
    let dir = std::env::temp_dir().join("rsomics-peak-count-compat");
    let _ = std::fs::create_dir_all(&dir);
    let bam = dir.join("in.bam");
    let peaks = golden("peaks.bed");
    // bedtools multicov needs a coordinate-sorted + indexed BAM; ours scans, so
    // sort once for a fair comparison (sorting doesn't change overlap counts).
    assert!(
        Command::new("samtools")
            .args(["sort", "-o"])
            .arg(&bam)
            .arg(golden("small.bam"))
            .status()
            .unwrap()
            .success()
    );
    assert!(
        Command::new("samtools")
            .arg("index")
            .arg(&bam)
            .status()
            .unwrap()
            .success()
    );

    let ours_out = String::from_utf8(
        Command::new(ours())
            .arg(&bam)
            .args(["-b", &peaks])
            .output()
            .unwrap()
            .stdout,
    )
    .unwrap();
    let bt_out = String::from_utf8(
        Command::new("bedtools")
            .args(["multicov", "-bed", &peaks, "-bams"])
            .arg(&bam)
            .output()
            .unwrap()
            .stdout,
    )
    .unwrap();

    // both put the count in the last tab-separated column
    let last_col = |s: &str| -> Vec<String> {
        s.lines()
            .filter(|l| !l.is_empty())
            .map(|l| l.rsplit('\t').next().unwrap().to_owned())
            .collect()
    };
    assert_eq!(
        last_col(&ours_out),
        last_col(&bt_out),
        "per-peak counts must match bedtools multicov"
    );
}
