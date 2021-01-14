import time

def GET_CUR_TIME(s=None):
	return "[{}]".format(time.asctime()) + '\t' + s if s else "[{}]".format(time.asctime()) 

def ROUND_DOWN(a, n):
	return a - a % n

def FLATTEN(x): return [i for s in x for i in s]

def ROUND_UP(a, n):
	return ROUND_DOWN(a + n - 1, n)

def PRINT_LOGO():
	print("""
 _____  __   _____  ___    __  ___      _ _           
/__   \\/__\\  \\_   \\/ __\\  /__\\/ __\\__ _| | | ___ _ __ 
  / /\\/ \\//   / /\\/__\\// /_\\ / /  / _` | | |/ _ \\ '__|
 / / / _  \\/\\/ /_/ \\/  \\//__/ /__| (_| | | |  __/ |   
 \\/  \\/ \\_/\\____/\\_____/\\__/\\____/\\__,_|_|_|\\___|_|                                               
	""")