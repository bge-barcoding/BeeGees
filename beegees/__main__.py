"""BeeGees command-line interface."""
import argparse
import shutil
import subprocess
import sys
import threading
from pathlib import Path

from beegees import __version__
from beegees.utils.configs import get_package_dir
from beegees.utils.snakemake_args import build_snakemake_cmd


def cmd_init(args):
    """Copy template config/ and profiles/ to the current working directory."""
    pkg = get_package_dir()
    cwd = Path.cwd()
    for d in ("config", "profiles"):
        dest = cwd / d
        if dest.exists() and not args.force:
            print(f"ERROR: {dest} already exists. Use --force to overwrite.", file=sys.stderr)
            return 1
        shutil.copytree(pkg / d, dest, dirs_exist_ok=args.force)
    samples_dest = cwd / "samples_template.csv"
    if not samples_dest.exists() or args.force:
        shutil.copy2(pkg / "config" / "samples_template.csv", samples_dest)
    print(f"Initialised BeeGees config in {cwd}")
    print(f"  1. Edit config/config.yaml")
    print(f"  2. Fill in samples_template.csv and save as samples.csv")
    print(f"  3. Run: beegees run --config config/config.yaml")
    return 0


def _reduce_thread_stack_size():
    """Reduce stack size for new Python threads to avoid virtual-memory exhaustion on HPC login nodes.

    Uses threading.stack_size() rather than resource.setrlimit(RLIMIT_STACK) so that
    the limit applies only to Python threads spawned within this process. Unlike
    setrlimit, threading.stack_size() is NOT inherited by child processes or SLURM
    jobs, preventing stack-overflow segfaults (exit 139) in tools like MitoGeneExtractor
    that rely on deep C/C++ call stacks on compute nodes.
    """
    try:
        threading.stack_size(2 * 1024 * 1024)  # 2 MB per Python thread
    except (threading.ThreadError, ValueError):
        pass


def cmd_run(args):
    """Run the BeeGees pipeline."""
    _reduce_thread_stack_size()
    configfile = Path(args.config)
    if not configfile.exists():
        print(f"ERROR: config file not found: {configfile}", file=sys.stderr)
        print(f"  Tip: run 'beegees init' to create a template config in the current directory,", file=sys.stderr)
        print(f"       or supply the path to an existing config with --config /path/to/config.yaml", file=sys.stderr)
        return 1

    log_file = Path(args.log_file) if args.log_file else None

    if not args.dryrun:
        unlock_cmd = build_snakemake_cmd(
            configfile=configfile,
            cores=args.cores,
            profile=args.profile,
            dryrun=False,
            unlock=True,
            extra_args=[],
        )
        result = subprocess.run(unlock_cmd, capture_output=True, text=True)
        if "Unlocked" in result.stdout or "Unlocked" in result.stderr:
            print("Working directory unlocked.")

    run_cmd = build_snakemake_cmd(
        configfile=configfile,
        cores=args.cores,
        profile=args.profile,
        dryrun=args.dryrun,
        unlock=False,
        extra_args=args.snakemake_args or [],
        log_file=log_file,
    )
    return subprocess.run(run_cmd).returncode


def build_parser():
    parser = argparse.ArgumentParser(
        prog="beegees",
        description="BeeGees — DNA barcoding from genome skims.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    # beegees init — optional convenience for CLI users
    p_init = sub.add_parser(
        "init",
        help="Copy template config/ and profiles/ to the current directory.",
    )
    p_init.add_argument(
        "--force", action="store_true",
        help="Overwrite existing config/ and profiles/ directories.",
    )
    p_init.set_defaults(func=cmd_init)

    # beegees run — Galaxy-compatible: --config is required
    p_run = sub.add_parser("run", help="Run the BeeGees pipeline.")
    p_run.add_argument(
        "--config", required=True, metavar="PATH",
        help="Path to config.yaml. Run 'beegees init' to generate a template.",
    )
    p_run.add_argument(
        "--cores", type=int, default=None, metavar="N",
        help="Number of CPU cores (overrides the profile setting).",
    )
    p_run.add_argument(
        "--profile", default=None, metavar="NAME_OR_PATH",
        help=(
            '"local" (default) or "slurm" use bundled profiles; '
            "or supply a path to a custom profile directory."
        ),
    )
    p_run.add_argument(
        "--dryrun", "-n", action="store_true",
        help="Preview jobs without executing.",
    )
    p_run.add_argument(
        "--log-file", default=None, metavar="PATH",
        help="Write the Snakemake log to this file (in addition to the terminal).",
    )
    p_run.add_argument(
        "snakemake_args", nargs=argparse.REMAINDER,
        help="Extra flags passed directly to snakemake (place after --).",
    )
    p_run.set_defaults(func=cmd_run)

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()
