import nox


nox.options.default_venv_backend = "uv"

@nox.session(python="3.13")
def test(session):
    session.install(".[test]")
    session.install("pytest-cov")
    session.run("uv", "pip", "list")
    session.run("pytest", "--durations=50", "tests", *session.posargs)



# 1) List your modules by hand (drop the "src/" and ".py"):
ALL_SCRIPTS = [
    "string_gsea.string_gsea_builder",
    "string_gsea.string_gsea_results",
    "string_gsea.gsea_result_processor",
]

@nox.session(name="run-internal-scripts")
def run_internal_scripts(session):
    # Install in editable mode
    session.install("-e", ".")

    # Determine which modules to run
    to_run = session.posargs or ALL_SCRIPTS
    if not to_run:
        session.error("No scripts to run!")

    # Record results here
    results = {}

    for module in to_run:
        session.log(f"▶ Running python -m {module}")
        try:
            # This will raise if the script returns non-zero
            session.run("python", "-m", module)
            results[module] = True
        except Exception:
            # Catch the failure, record it, and keep going
            results[module] = False

    # Print a summary
    session.log("")
    session.log("▶ Script run summary:")
    for module, ok in results.items():
        mark = "✓" if ok else "✗"
        session.log(f"  {mark} {module}")

    # If any failed, signal overall failure
    if not all(results.values()):
        session.error("One or more scripts failed.")

