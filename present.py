import matplotlib.pyplot as plt

# Needs to become a class
data_dict: dict[tuple | tuple[str, int]] = {}
data_dict[(5, 5)] = ("P", 0)
data_dict[(5, 6)] = ("H", 1)
data_dict[(4, 6)] = ("P", 2)
data_dict[(4, 7)] = ("P", 3)
data_dict[(5, 7)] = ("H", 4)
data_dict[(6, 7)] = ("P", 5)
data_dict[(6, 6)] = ("H", 6)
data_dict[(6, 5)] = ("H", 7)

# Needs to become another class
def show_points(data: dict[tuple | tuple[str, int]]) -> None:
    # Filter the points out of the data
    # Polair
    x_p = [item[0][0] for item in data.items() if item[1][0] == "P"]
    y_p = [item[0][1] for item in data.items() if item[1][0] == "P"]

    # Hydrofobe
    x_h = [item[0][0] for item in data.items() if item[1][0] == "H"]
    y_h = [item[0][1] for item in data.items() if item[1][0] == "H"]

    # Covalent bonds
    x_l = [item[0][0] for item in data.items()]
    y_l = [item[0][1] for item in data.items()]

    # Borders for plot
    x_min, x_max, y_min, y_max = min(x_l), max(x_l), min(y_l), max(y_l)

    # Make the plot with dots and line
    plt.plot(x_l, y_l, c = "black", alpha= 0.8, linewidth= 5)
    plt.plot(x_p, y_p, "bo", markersize= 20)
    plt.plot(x_h, y_h, "ro", markersize= 20)
    plt.xlim(x_min - 1, x_max + 1)
    plt.ylim(y_min - 1, y_max + 1)
    plt.legend(["Bond", "Polair", "Hydrofobe"])
    plt.show()

def options(y: int, x: int, data: dict[tuple | tuple[str, int]]) -> list[tuple[int]]:
    result = []

    # Define the bounds
    # x-1, y
    # x  , y+1 
    # x+1, y
    # x  , y-1
    

    for i in range(-1, 3):
        print((i % 2) * i // abs(i), 1 - abs(i))
        # if data[()]
        # for j in range(2):
    #         diff = abs(y - i ) + abs(x - j)
    #         if diff == 1:
    #             result.append((i, j))
    return result

def calculate_stability(data: list[list[str | tuple[str, int]]]) -> int:
    filtered_list = list(filter(lambda x: x[0] == "H", filtered_list))
    filtered_dict = {}
    score = 0
    for value in filtered_list:
        filtered_dict[f"{value[2], value[3]}"] = value[0]
    for value in filtered_list:
        for option in options(value[2], value[3]):
            if f"{option}" in filtered_dict and filtered_dict[f"{option}"] == "H":
                print("ja")

options(3, 5, data_dict)