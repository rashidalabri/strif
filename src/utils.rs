use std::path::PathBuf;

pub struct AlignmentScoreParams {
    pub match_score: i32,
    pub mismatch_penalty: i32,
    pub gap_open_penalty: i32,
    pub gap_extend_penalty: i32,
}

pub fn get_default_out_path(input: &PathBuf, suffix: &str, ext: &str) -> PathBuf {
    let mut out_path: PathBuf = input.clone();
    let mut file_prefix = input.file_stem().unwrap().to_str().unwrap();

    // extract text before first period
    if let Some(period_idx) = file_prefix.find('.') {
        file_prefix = &file_prefix[..period_idx];
    }

    out_path.set_file_name(format!("{}.{}.{}", file_prefix, suffix, ext));
    out_path
}
