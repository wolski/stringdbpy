import nox


nox.options.default_venv_backend = "uv"

@nox.session(python="3.13")
def test(session):
    session.install(".[test]")
    session.run("uv", "pip", "list")
    session.run("pytest", "--durations=50", "tests", *session.posargs)
