use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_peak_count(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-peak-count");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let bam = manifest.join("tests/golden/small.bam");
    let bed = manifest.join("tests/golden/peaks.bed");

    c.bench_function("rsomics-peak-count golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args([bam.to_str().unwrap(), "--bed", bed.to_str().unwrap()])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_peak_count);
criterion_main!(benches);
