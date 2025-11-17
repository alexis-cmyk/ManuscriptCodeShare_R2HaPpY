# 99_update_and_snapshot.R
# Run this ONLY on your "good" machine when you want to refresh all package versions.

# Make sure we're in the project and renv is active
if (Sys.getenv("RENV_PROJECT") == "") {
  stop("renv is not active – open the project and try again.")
}

# 1) Ignore base + recommended packages (Matrix, boot, DBI, etc.)
base_recommended <- rownames(installed.packages(priority = c("base", "recommended")))

current_ignored <- renv::settings$ignored.packages()
renv::settings$ignored.packages(
  sort(unique(c(current_ignored, base_recommended)))
)

message("Ignoring base/recommended packages:\n",
        paste(renv::settings$ignored.packages(), collapse = ", "))

# 2) Update ALL non-ignored packages to the latest compatible versions
#    (CRAN + Bioc, using your repos from .Rprofile)
message("Updating all non-ignored packages via renv::update() ...")
renv::update(prompt = FALSE)

# 3) Snapshot the result into renv.lock
message("Snapshotting updated environment to renv.lock ...")
renv::snapshot(type = "all")

message("Done. Commit the updated renv.lock (and renv/settings.dcf if modified).")