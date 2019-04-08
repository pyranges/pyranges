def test_tssify(f1):

    result = f1.tssify(slack=5)

    print(f1)
    print(result)

    assert list(result.Start) == [0, 3, 2]
    assert list(result.End) == [9, 14, 13]


def test_tesify(f2):

    print(f2)

    result = f2.tesify(slack=500)

    print(result)

    assert list(result.Start) == [0, 0]
    assert list(result.End) == [503, 508]
