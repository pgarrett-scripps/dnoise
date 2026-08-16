//! dnoise desktop GUI — a point-and-click front end over the `dnoise` library for
//! users who would rather not touch the command line. Drag `.d` folders in, pick an
//! acquisition preset and an output location, and Run.

// On Windows, don't pop a console window behind the GUI in release builds.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod app;
mod settings;
mod worker;

use app::DnoiseApp;

fn main() -> eframe::Result<()> {
    let native_options = eframe::NativeOptions {
        viewport: eframe::egui::ViewportBuilder::default()
            .with_inner_size([740.0, 680.0])
            .with_min_inner_size([560.0, 480.0])
            .with_title("dnoise"),
        ..Default::default()
    };
    eframe::run_native(
        "dnoise",
        native_options,
        Box::new(|_cc| {
            let mut app = DnoiseApp::default();
            // Directories passed on the command line (Explorer context menu,
            // drag-onto-exe) are queued as inputs. add_path validates that each
            // is a .d folder and logs a skip note otherwise.
            for arg in std::env::args_os().skip(1) {
                let p = std::path::PathBuf::from(arg);
                if p.is_dir() {
                    app.add_path(p);
                }
            }
            Ok(Box::new(app))
        }),
    )
}
