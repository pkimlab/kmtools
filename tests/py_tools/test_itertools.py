import itertools

from kmtools import py_tools

# === iter_forever ===


def test_iter_forever_1():

    def my_generator():
        yield from range(3)

    g = py_tools.iter_forever(my_generator)

    list(itertools.islice(g, 7)) == [0, 1, 2, 0, 1, 2, 0]


def test_iter_forever_2():

    def my_generator():
        yield from range(3)

    g = py_tools.iter_forever(my_generator)

    outputs = []
    for i in range(1_000):
        outputs.append(next(g))

    assert outputs == list(itertools.islice(itertools.cycle(my_generator()), 1_000))


def test_iter_forever_3():

    def my_generator():
        i = 0
        for j in range(3):
            i = yield f"{i} {j}"

    g = py_tools.iter_forever(my_generator)
    g.send(None)

    for i in itertools.islice(itertools.cycle(range(3)), 1, 1_000):
        value = g.send(i)
        assert value == f"{i} {i}"
