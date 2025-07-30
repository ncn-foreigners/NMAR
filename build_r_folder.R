# build_r_folder.R

# --- Configuration ---
source_directory <- "src"  # The folder to copy files from (your source code)
target_directory <- "R"    # The destination folder for R package files

# --- Helper function to generate the target filename ---
# Takes the full path to a source file and the base source directory.
# Returns a new filename like "models_linear.R" from "src/models/linear.R".
generate_target_filename <- function(file_path, base_source_dir) {
  # Remove the base source directory from the beginning of the path
  relative_path <- sub(paste0("^", base_source_dir, "/?"), "", file_path)

  # Replace "/" with "_" in the path and remove the ".R" extension.
  # Then, add ".R" back at the end.
  new_name <- gsub("/", "_", relative_path)
  new_name <- gsub("\\.R$", "", new_name) # Remove .R from the end
  new_name <- paste0(new_name, ".R")      # Add .R back

  return(new_name)
}

# --- Main Script Logic ---

message("--- Preparing the R/ folder ---")

# 1. Check if the source directory exists
if (!dir.exists(source_directory)) {
  stop(paste0("Error: Source directory '", source_directory, "' does not exist."))
}

# 2. Create the target directory (R/) if it doesn't exist
if (!dir.exists(target_directory)) {
  dir.create(target_directory)
  message(paste0("Created target folder: ", target_directory, "/"))
}

# 3. Remove all existing .R files from the target directory to avoid duplicates or outdated files
existing_r_files <- list.files(target_directory, pattern = "\\.R$", full.names = TRUE)
if (length(existing_r_files) > 0) {
  file.remove(existing_r_files)
  message(paste0("Removed ", length(existing_r_files), " existing .R files from ", target_directory, "/"))
}

# 4. Find all .R files in the source directory and its subdirectories
all_source_r_files <- list.files(
  source_directory,
  pattern = "\\.R$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(all_source_r_files) == 0) {
  message(paste0("No .R files found in '", source_directory, "/'."))
} else {
  message(paste0("Found ", length(all_source_r_files), " .R files to copy."))
  # 5. Copy files to the target directory with renamed paths
  for (src_file_path in all_source_r_files) {
    target_filename <- generate_target_filename(src_file_path, source_directory)
    target_file_path <- file.path(target_directory, target_filename)

    file.copy(src_file_path, target_file_path, overwrite = TRUE)
    message(paste0("Copied: ", src_file_path, " -> ", target_file_path))
  }
}

message("--- Finished preparing the R/ folder ---")
