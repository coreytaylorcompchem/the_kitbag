import random

def expand_seeds(num_runs, base_seed=None, strategy="increment"):
    if num_runs < 1:
        raise ValueError("num_runs must be >= 1")

    if strategy == "fixed":
        if base_seed is None:
            raise ValueError("base_seed required for fixed strategy")
        return [base_seed] * num_runs

    if strategy == "increment":
        if base_seed is None:
            base_seed = random.randint(0, 2**31 - 1)
        return [base_seed + i for i in range(num_runs)]

    if strategy == "random":
        return [random.randint(0, 2**31 - 1) for _ in range(num_runs)]

    raise ValueError(f"Unknown seed strategy: {strategy}")
