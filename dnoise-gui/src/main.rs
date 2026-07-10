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
        Box::new(|_cc| Ok(Box::new(DnoiseApp::default()))),
    )
}
