# -*- coding: utf-8 -*-
#!/usr/bin/env python

import time
from TRIBECaller._version import __version__


def GET_CUR_TIME(s=None):
    return "[{}]".format(time.asctime()) + '\t' + s if s else "[{}]".format(time.asctime()) 

def ROUND_DOWN(a, n):
    return a - a % n

def FLATTEN(x): return [i for s in x for i in s]

def CONDENSE(arr, bin_size):
    res = []
    temp = 0
    for i in range(len(arr)):
        if i % bin_size != bin_size-1:
            temp += arr[i]
        else:
            res.append(temp+arr[i])
            temp=0
    return res

def PICK(arr, bin_size):
    res = []
    for i in range(len(arr)):
        if i % bin_size == bin_size-1:
            res.append(arr[i])
    return res

def ROUND_UP(a, n):
    return ROUND_DOWN(a + n - 1, n)

def PRINT_PRELOGUE():
    print("")
    print("Welcome to")

def PRINT_LOGO():
    print(GET_RED(""" _____  __   _____  ___    __  ___      _ _   """))
    print(GET_BLUE("""/__   \\/__\\  \\_   \\/ __\\  /__\\/ __\\__ _| | | ___ _ __ """))
    print(GET_YELLOW("""  / /\\/ \\//   / /\\/__\\// /_\\ / /  / _` | | |/ _ \\ '__|"""))
    print(GET_GREEN(""" / / / _  \\/\\/ /_/ \\/  \\//__/ /__| (_| | | |  __/ |   """))
    print(GET_PURPLE(""" \\/  \\/ \\_/\\____/\\_____/\\__/\\____/\\__,_|_|_|\\___|_| """))
    print("")                                       

def PRINT_INFO():
    print("Version: " + __version__)

def GET_RED(s):
    return "\x1b[{};{};{}m".format(BackgroundLighness.LEVEL2, FontColor.RED, BackgroundColor.NULL) + s + "\x1b[0m"

def GET_GREEN(s):
    return "\x1b[{};{};{}m".format(BackgroundLighness.LEVEL2, FontColor.GREEN, BackgroundColor.NULL) + s + "\x1b[0m"

def GET_YELLOW(s):
    return "\x1b[{};{};{}m".format(BackgroundLighness.LEVEL2, FontColor.YELLOW, BackgroundColor.NULL) + s + "\x1b[0m"

def GET_BLUE(s):
    return "\x1b[{};{};{}m".format(BackgroundLighness.LEVEL2, FontColor.BLUE, BackgroundColor.NULL) + s + "\x1b[0m"

def GET_PURPLE(s):
    return "\x1b[{};{};{}m".format(BackgroundLighness.LEVEL2, FontColor.PURPLE, BackgroundColor.NULL) + s + "\x1b[0m"

def GET_JASPER(s):
    return "\x1b[{};{};{}m".format(BackgroundLighness.LEVEL2, FontColor.JASPER, BackgroundColor.NULL) + s + "\x1b[0m"


class BackgroundLighness:
    LEVEL1 = 0
    LEVEL2 = 1
    LEVEL3 = 2
    LEVEL4 = 3
    LEVEL5 = 4
    LEVEL6 = 5
    LEVEL7 = 6
    REVERSE= 7

class FontColor:
    BLACK = 30
    RED = 31
    GREEN = 32
    YELLOW = 33
    BLUE = 34
    PURPLE = 35
    JASPER = 36
    WHITE = 37

class BackgroundColor:
    BLACK = 40
    RED = 41
    GREEN = 42
    YELLOW = 43
    BLUE = 44
    PURPLE = 45
    JASPER = 46
    WHITE = 47
    NULL = 48

NUC = {0:('A','#00B500'),1:('T','#FF005C'),2:('C','#0083FF'),3:('G','#FFA700')}