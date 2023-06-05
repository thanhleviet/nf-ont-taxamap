/**
 * Scans a folder for files and folders matching specific patterns and returns a list of tuples.
 * Each tuple contains two elements: a string and a boolean value.
 * The string represents a wildcard pattern that can be used to match files or folders.
 * The boolean value is true if the string represents a folder pattern, false if it represents a file pattern.
 * If no files or folders matching the specified patterns are found, the function returns a tuple with the string "none" and the value false.
 *
 * @param folderPath The path to the folder to scan.
 * @return A list of tuples representing the files and folders matching the specified patterns.
 */
def scanFolder(folderPath) {
    def results = []
    def folder = new File(folderPath)
    if (folder.exists() && folder.isDirectory()) {
        folder.eachFile { file ->
            if (file.isDirectory() && file.name.matches(/barcode\d+/)) {
                // If the file is a directory and its name matches the pattern "barcode\d+",
                // add a tuple to the results list with the string "barcode*" and the value true.
                results << ["dir","barcode*"]
            } else if (file.isFile() && (file.name.matches(/barcode.*\.fastq\.gz/) || file.name.matches(/barcode.*\.fq\.gz/))) {
                // If the file is a regular file and its name matches the pattern "barcode.*\.fastq\.gz" or "barcode.*\.fq\.gz",
                // add a tuple to the results list with the string "barcode*.fastq.gz" or "barcode*.fq.gz" and the value false.
                if (file.name.matches(/barcode.*\.fastq\.gz/)) {
                    results << ["file","barcode*.fastq.gz"]
                } else {
                    results << ["file","barcode*.fq.gz"]
                }
            }
        }
    }
    if (results.isEmpty()) {
        // If no files or folders matching the specified patterns are found,
        // add a tuple to the results list with the string "none" and the value false.
        results << ["none", false]
    }
    return results
}