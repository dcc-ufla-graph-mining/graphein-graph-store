import pickle
from memory_profiler import profile

@profile
def main():
    with open('v1.pkl', 'rb') as f:
        v1 = pickle.load(f)


if __name__ == "__main__":
    main()
