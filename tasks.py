from invoke import task


@task
def tidy(c, docs=False):
    c.run("isort . --skip .conda")
    c.run("black . --exclude .conda")
    c.run("flake8 . --exclude .conda")