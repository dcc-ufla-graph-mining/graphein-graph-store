import pickle
from memory_profiler import profile

@profile
def main():
    print(pickle.__name__)

main()