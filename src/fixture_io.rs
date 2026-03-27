use flate2::read::GzDecoder;
use std::ffi::OsString;
use std::fs;
use std::fs::File;
use std::io;
use std::io::Read;
use std::path::{Path, PathBuf};

fn gz_candidate(path: &Path) -> PathBuf {
    let mut candidate: OsString = path.as_os_str().to_os_string();
    candidate.push(".gz");
    PathBuf::from(candidate)
}

pub fn resolve_fixture_path(path: impl AsRef<Path>) -> Option<PathBuf> {
    let path = path.as_ref();
    if path.exists() {
        return Some(path.to_path_buf());
    }

    let gz_path = gz_candidate(path);
    if gz_path.exists() {
        return Some(gz_path);
    }

    None
}

pub fn fixture_exists(path: impl AsRef<Path>) -> bool {
    resolve_fixture_path(path).is_some()
}

pub fn read_text_fixture(path: impl AsRef<Path>) -> io::Result<String> {
    let requested_path = path.as_ref();
    let resolved_path = resolve_fixture_path(requested_path).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "fixture not found: {} or {}.gz",
                requested_path.display(),
                requested_path.display()
            ),
        )
    })?;

    if resolved_path.extension().and_then(|ext| ext.to_str()) == Some("gz") {
        let file = File::open(&resolved_path)?;
        let mut decoder = GzDecoder::new(file);
        let mut content = String::new();
        decoder.read_to_string(&mut content)?;
        return Ok(content);
    }

    fs::read_to_string(resolved_path)
}
