#!/usr/bin/env python3

import pytest
import pyranges as pr

def assert_equal_length_before_after(gr1, gr2):

    print("in test")
    l1 = len(gr1)
    l2 = len(gr2)
    print(pr.concat([gr1, gr2]))
    c = len(pr.concat([gr1, gr2]))

    assert l1 + l2 == c

def test_concat_stranded_unstranded(f1, f2):

    assert_equal_length_before_after(f1, f2)

def test_concat_unstranded_unstranded(f1, f2):

    assert_equal_length_before_after(f1.unstrand(), f2.unstrand())

def test_concat_stranded_unstranded(f1, f2):
    assert_equal_length_before_after(f1, f2.unstrand())

def test_concat_unstranded_stranded(f1, f2):
    assert_equal_length_before_after(f1.unstrand(), f2)
