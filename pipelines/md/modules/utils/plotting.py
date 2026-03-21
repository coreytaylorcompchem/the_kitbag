import matplotlib.pyplot as plt

def save_plot(x, y, xlabel, ylabel, title, filepath, labels=None):
    plt.figure()

    if isinstance(y, dict):  # multiple lines
        for key, values in y.items():
            plt.plot(x, values, label=key)
        if labels:
            plt.legend()
    else:
        plt.plot(x, y)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.savefig(filepath)
    plt.close()