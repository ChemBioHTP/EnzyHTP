"""Testing the enzy_htp.core.EnvironmentManager() object.

Author: Chris Jurich <chris.jurich@vanderbilt.edu>
Date: 2022-03-19
"""
import os
import pytest

from enzy_htp.core import EnvironmentManager


def test_check_environment():
    """Ensuring that the manager can check a dummy environment."""
    os.environ["___VAR1___"] = "VALUE1"
    os.environ["___VAR2___"] = "VALUE2"
    os.environ["___VAR3___"] = "VALUE3"
    em = EnvironmentManager()
    em.add_env_var("___VAR1___")
    em.add_env_var("___VAR2___")
    em.add_env_var("___VAR3___")
    em.add_executable("ls")
    em.add_executable("less")
    em.add_executable("cat")
    em.check_environment()
    assert not em.is_missing()

    em = EnvironmentManager(
        env_vars=["___VAR1___", "___VAR2___", "___VAR3___"],
        executables=["ls", "less", "cat"],
    )
    em.check_environment()
    assert not em.is_missing()
    os.unsetenv("___VAR1___")
    os.unsetenv("___VAR2___")
    os.unsetenv("___VAR3___")


def test_missing_env_var():
    """Ensuring that the manager effectivetly detects a missing environment variable."""
    em = EnvironmentManager()
    em.add_env_var("___VAR4___")
    em.check_environment()
    assert em.is_missing()
    em.reset()
    assert not em.is_missing()


def test_missing_executable():
    """Ensuring that the manager effectively detects a missing exceutable."""
    em = EnvironmentManager()
    em.add_executable("does-not-exist")
    em.check_environment()
    assert em.is_missing()
    em.reset()
    assert not em.is_missing()


def test_executable___gettattr__():
    """Ensuring that variable accession works for acquired variables."""
    em = EnvironmentManager(executables=["ls"])
    with pytest.raises(SystemExit) as exe:
        em.ls
    assert exe.type == SystemExit
    assert exe.value.code == 1

    em.check_environment()
    assert len(em.ls)

    with pytest.raises(KeyError) as exe:
        em.dne

    assert str(exe).find("dne") != -1


def test_run_command():
    """Ensuring that EnvironmentManager.run_command() works and raises errors when it is supposed to."""
    em = EnvironmentManager()
    contents_str = em.run_command("echo", ["hello world"])
    assert contents_str == ['"hello world"']
    
    contents_file_name = em.run_command("echo", ["test.txt"])
    assert contents_file_name == ['"test.txt"']

    contents_special_chr = em.run_command("echo", [">", "&&"])
    assert contents_special_chr == ['> &&']

    contents_all_kinds = em.run_command("echo", ["input.inp", "I am Enzyme", ">", "&&" , ">>>!@#$%^&*()"])
    assert contents_all_kinds == ['"input.inp" "I am Enzyme" > && >>>!@#$%^&*()']

    with pytest.raises(SystemExit) as exe:
        em.run_command("dne", [""])

    assert exe.type == SystemExit
    assert exe.value.code == 1

    em.check_environment()

    with pytest.raises(SystemExit) as exe:
        em.run_command("dne", [""])

    assert exe.type == SystemExit
    assert exe.value.code == 1
