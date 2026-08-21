//! Support for the gzip compressed bincode databases that crates building on top of this one embed
//! in the binary at compile time.

use std::io::Read;

/// Decompress a gzip compressed static database.
///
/// The databases that ship inside the binary are stored compressed because the uncompressed bincode
/// is dominated by long range redundancy and dwarfs the actual code (GNOme alone is 25 MB
/// uncompressed against 8 MB compressed). Decompression happens once, behind the `LazyLock` that
/// guards the data, and adds roughly a fifth to the bincode decode that has to happen regardless.
///
/// This is deliberately kept separate from the runtime cache, which stays uncompressed, so that
/// caches written by older versions stay readable.
///
/// # Panics
/// If the data is not valid gzip. This can only happen if the database that was embedded at compile
/// time is corrupt.
#[must_use]
pub fn decompress_static_data(compressed: &[u8]) -> Vec<u8> {
    let mut decompressed = Vec::new();
    flate2::read::GzDecoder::new(compressed)
        .read_to_end(&mut decompressed)
        .expect("An embedded static database is not valid gzip data");
    decompressed
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::decompress_static_data;

    /// The databases are written by `CVIndex::save_to_static_data_at` with a
    /// `flate2::write::GzEncoder` and read back with [`decompress_static_data`]. Nothing else pins
    /// those two together, so a change to either side that breaks the pairing would only show up
    /// the next time the databases are regenerated.
    #[test]
    fn encoder_and_decoder_are_inverse() {
        let data: Vec<u8> = (0..10_000_u32).map(|i| (i % 251) as u8).collect();

        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::best());
        encoder.write_all(&data).unwrap();
        let compressed = encoder.finish().unwrap();

        assert!(
            compressed.len() < data.len(),
            "the test data should compress"
        );
        assert_eq!(decompress_static_data(&compressed), data);
    }
}
