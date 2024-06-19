import nox


# An example nox task definition that runs on many supported Python versions:
@nox.session(
    python=["3.7", "3.8", "3.9", "3.10", "3.11", "3.12"]
)
def test(session):
    session.install(".")

    session.run("python", "-m", "unittest")


@nox.session
@nox.parametrize("python,numpy", [("3.7", "1.16"), ("3.9", "1.20"), ("3.10", "2.0")])
def numpy(session, numpy):
    session.install(f"numpy=={numpy}")
    session.install(".")

    session.run("python", "-m", "unittest")


@nox.session
@nox.parametrize("python,plotly", [("3.7", "5.12"), ("3.9", "5.17"),  ("3.10", "5.22")])  # 5.0 does not work with newer (>1.23) numpy due to np.bool deprecation
def plotly(session, plotly):
    session.install(f"plotly=={plotly}")
    session.install(".")

    session.run("python", "-m", "unittest")


@nox.session
@nox.parametrize("python,matplotlib", [("3.7", "3.1"), ("3.8", "3.3"), ("3.10", "3.9")])
def matplotlib(session, matplotlib):
    session.install(f"matplotlib=={matplotlib}")
    session.install(".")

    session.run("python", "-m", "unittest")
