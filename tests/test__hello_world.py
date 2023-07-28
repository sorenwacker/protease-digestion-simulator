from template import hello_world


def test__hello_world():
  expected = 'Hello world!'
  actual = hello_world()
  assert actual == expected, f"hello_world() returned '{actual}' instead of '{expected}'."
  
