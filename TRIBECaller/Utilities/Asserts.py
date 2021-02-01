# -*- coding: utf-8 -*-
#!/usr/bin/env python
# Author: Ziwei Xue
#
# Custom Assertion Functions

from TRIBECaller.Utilities.utils import *

def assert_eq(a,b):
	if a != b:
		raise AssertionError(GET_RED("Assertion Equality Failed!"))

def assert_geq():
	if a < b:
		raise AssertionError(GET_RED("Assertion Greater or Equal Failed!"))

def assert_leq():
	if a > b:
		raise AssertionError(GET_RED("Assertion Less or Equal Failed!"))

def assert_in():
	if a not in b:
		raise AssertionError(GET_RED("Assertion in Failed!"))