def test_five_end(f1):

    result = f1.five_end(slack=5)

    print(f1)
    print(result)

    assert list(result.Start) == [0, 3, 2]
    assert list(result.End) == [9, 14, 13]


def test_three_end(f2):

    print(f2)

    result = f2.three_end(slack=500)

    print(result)

    assert list(result.Start) == [0, 0]
    assert list(result.End) == [503, 508]
