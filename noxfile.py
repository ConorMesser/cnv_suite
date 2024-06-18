import nox


# An example nox task definition that runs on many supported Python versions:
@nox.session(
    python=["3.7", "3.8", "3.9", "3.10", "3.11", "3.12"]
)
def test(session):
    session.install(".")

    session.run("python", "-m", "unittest")


@nox.session()
@nox.parametrize("numpy", ["1.14", "2.0"])
def numpy(session, numpy):
    session.install(f"numpy=={numpy}")
    session.install(".")

    session.run("python", "-m", "unittest")


@nox.session()
@nox.parametrize("plotly", ["5.0", "5.22"])
def plotly(session, plotly):
    session.install(f"plotly=={plotly}")
    session.install(".")

    session.run("python", "-m", "unittest")


@nox.session()
@nox.parametrize("matplotlib", ["3.0", "3.9"])
def matplotlib(session, matplotlib):
    session.install(f"matplotlib=={matplotlib}")
    session.install(".")

    session.run("python", "-m", "unittest")
