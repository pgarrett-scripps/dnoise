//! Parity test for the streaming API ([`dnoise::RunContext`]).
//!
//! The file writer ([`dnoise::denoise`]) and the streaming context share the same
//! per-frame core ([`dnoise::process_frame_decoded`]); the only thing that could
//! drift is the per-run gate *wiring* each assembles. This test denoises a real
//! `.d` with the writer, reads the survivors back out of the written `.d` (via the
//! same timsrust decode dnoise trusts), and asserts they match — frame by frame —
//! what `RunContext::process` returns on the same input.
//!
//! It needs a real `.d`, so it is gated on the `DNOISE_PARITY_DOTD` env var and
//! is a no-op when that is unset (keeps CI green without a fixture):
//!
//! ```sh
//! DNOISE_PARITY_DOTD=/path/to/example.d cargo test --test streaming -- --nocapture
//! ```

use std::path::PathBuf;

use dnoise::frame::FlatFrame;
use dnoise::{FilterParams, RunContext, Stages, denoise};
use timsrust::readers::FrameReader;

/// Sort a survivor set into a canonical order for comparison (the encoder groups
/// by scan and sorts by TOF, so raw order is not meaningful).
fn canonical(mut pts: Vec<(u32, u32, u32)>) -> Vec<(u32, u32, u32)> {
    pts.sort_unstable();
    pts
}

fn flat_points(flat: &FlatFrame) -> Vec<(u32, u32, u32)> {
    (0..flat.len())
        .map(|j| (flat.scan[j], flat.tof[j], flat.intensity[j]))
        .collect()
}

// Ignored by default so a missing fixture shows up as an *ignored* test rather
// than a passing no-op. Run it with the fixture via:
//   DNOISE_PARITY_DOTD=/path/to/example.d cargo test --test streaming -- --ignored
#[test]
#[ignore = "requires DNOISE_PARITY_DOTD pointing at a real .d"]
fn streaming_matches_writer() {
    let input = match std::env::var("DNOISE_PARITY_DOTD") {
        Ok(p) => PathBuf::from(p),
        Err(_) => panic!(
            "DNOISE_PARITY_DOTD unset — set it to a real .d to run the writer/streaming parity test"
        ),
    };

    // Default denoising (vertical filter + halo on MS1). No crop, so the writer's
    // output frames line up 1:1 with RunContext frames.
    let params = FilterParams::default();
    let stages = Stages::default();

    let out = std::env::temp_dir().join("dnoise_parity_out.d");
    denoise(&input, &out, &params, &stages, true).expect("writer denoise");

    let ctx = RunContext::open(&input, &params, &stages).expect("open streaming context");
    let out_reader = FrameReader::new(&out).expect("open written .d");

    for i in 0..ctx.len() {
        let streamed = canonical(ctx.process(i).expect("stream frame").survivors);

        // timsrust cannot decode an absent (empty) payload; such a frame must have
        // produced no survivors on the streaming side too.
        match out_reader.get(i) {
            Ok(frame) => {
                let written = canonical(flat_points(&FlatFrame::from_frame(&frame)));
                assert_eq!(
                    written, streamed,
                    "frame {i}: writer survivors differ from RunContext survivors"
                );
            }
            Err(_) => assert!(
                streamed.is_empty(),
                "frame {i}: writer emitted empty but RunContext kept {} points",
                streamed.len()
            ),
        }
    }
}
