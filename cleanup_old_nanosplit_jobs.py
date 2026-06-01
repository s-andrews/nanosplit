#!/usr/bin/env python3
#
# Remove NanoSplit jobs older than 7 days.
#
# Linux cron install example:
#   cd /opt/nanosplit
#   chmod +x cleanup_old_jobs.py
#   sudo crontab -e
#
# Add this line to run the cleanup every day at 03:15:
#   15 3 * * * /opt/nanosplit/cleanup_old_jobs.py >> /var/log/nanosplit-cleanup.log 2>&1
#
# The script resolves the data directory using the same rule as app.py:
# the repository/application directory containing this script, plus "data".

import shutil
import time
from pathlib import Path


RETENTION_SECONDS = 7 * 24 * 60 * 60


def data_dir() -> Path:
    return Path(__file__).resolve().parent / "data"


def job_age_time(job_dir: Path) -> float:
    job_file = job_dir / "job.json"
    if job_file.exists():
        return job_file.stat().st_mtime
    return job_dir.stat().st_mtime


def main() -> int:
    root = data_dir()
    if not root.exists():
        print(f"Data directory does not exist: {root}")
        return 0
    if not root.is_dir():
        print(f"Data path is not a directory: {root}")
        return 1

    cutoff = time.time() - RETENTION_SECONDS
    removed = 0
    for job_dir in root.iterdir():
        if not job_dir.is_dir():
            continue
        if job_age_time(job_dir) >= cutoff:
            continue

        shutil.rmtree(job_dir)
        removed += 1
        print(f"Removed old job: {job_dir}")

    print(f"Removed {removed} old job(s) from {root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
