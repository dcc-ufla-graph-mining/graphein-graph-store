import pickle
from memory_profiler import profile

@profile
def main():
    with open('v2.pkl', 'rb') as f:
        v2 = pickle.load(f)



if __name__ == "__main__":
    main()