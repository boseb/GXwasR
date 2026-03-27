## Verify GCTA
verifyGCTA <- function() {
    os_type <- detect_os_type()
    # 1. Check GCTA_PATH environment variable
    env_path <- Sys.getenv("GCTA_PATH", unset = NA)
    if (!is.na(env_path)) {
        resolved_env_path <- normalizePath(env_path, mustWork = FALSE)
        if (file.exists(resolved_env_path)) {
            return(resolved_env_path)
        } else {
            message("GCTA_PATH is set but file does not exist: ", resolved_env_path)
        }
    }

    # 2. Check system path via Sys.which
    sys_path <- Sys.which("gcta64")
    # On Windows, check for `.exe` if missing
    if (os_type == "windows" && !grepl("\\.exe$", sys_path, ignore.case = TRUE)) {
        sys_path <- paste0(sys_path, ".exe")
    }
    # Normalize and check existence
    resolved_sys_path <- normalizePath(sys_path, mustWork = FALSE)
    if (nzchar(resolved_sys_path) && file.exists(resolved_sys_path)) {
        return(resolved_sys_path)
    }

    # 3. Failure message
    rlang::abort(
        message = rlang::format_error_bullets(c(
            "GCTA binary not found.",
            "x" = "Attempted to locate the 'GCTA_PATH' environment variable and 'gcta64' in system PATH.",
            "i" = "Ensure GCTA is available and executable.",
            "i" = "You can permanently set the path to 'gcta64' by running:",
            " " = "usethis::edit_r_environ()  # then add a line like: GCTA_PATH=/full/path/to/gcta64"
        )),
        class = "gcta_not_found"
    )
}

gcta <- function() {
    gcta_executable <- verifyGCTA()
    gcta_info <- suppressWarnings(system2(gcta_executable, stdout = TRUE))
    rlang::inform(
        message = paste("Using GCTA", stringr::str_remove(gcta_info[[3]], "\\* ")),
        .frequency = "regularly",
        .frequency_id = "gcta_check"
    )
    gcta_executable
}

## Added in 3.0
executeGCTA <- function(ResultDir, args) {
    tryCatch(
        {
            # Determine the destinations for standard output and standard error based on the operating system
            stdout_dest <- ifelse(.Platform$OS.type == "windows", "NUL", "/dev/null")
            stderr_dest <- stdout_dest # Redirecting stderr to the same null device

            file_path <- normalizePath(file.path(ResultDir, "sink_file.txt"), mustWork = FALSE)

            # Create an empty file (if it doesn't exist)
            if (!file.exists(file_path)) {
                file.create(file_path)
            }

            # Redirect output to the specified file
            sink(file_path)

            # Execute GCTA command

            invisible(sys::exec_wait(
                gcta(), # Path to the GCTA executable
                args = args, # Arguments for the GCTA command
                std_out = stdout_dest, # Standard output redirection
                std_err = stderr_dest # Standard error redirection
                # std_out = TRUE, # Standard output redirection
                # std_err = TRUE
            ))
        },
        error = function(e) {
            stop("An error occurred while executing GCTA: ", e$message)
        }
    )

    # Reset the sink to stop redirecting the output to the file
    if (sink.number() > 0) {
        sink() # This line resets the output redirection
    }
}