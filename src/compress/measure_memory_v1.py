import pickle
import time
from memory_profiler import profile

@profile
def main():
    with open('v1.pkl', 'rb') as f:
        v1 = pickle.load(f)
    
    time.sleep(3)


if __name__ == "__main__":
    main()
