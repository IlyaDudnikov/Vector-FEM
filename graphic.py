import matplotlib.pyplot as plt
import numpy as np

# Функция для чтения данных из файла
def read_csv(filename):
    return np.loadtxt(filename, delimiter=' ')

def vector_function(X, Y):
    U = np.cos(Y) / 10
    V = np.sin(X) / 10
    return U, V
# Чтение данных из файлов
points = read_csv('points.txt')
vectors = read_csv('vectors.txt')

# Проверка, что количество точек и векторов совпадает
assert points.shape[0] == vectors.shape[0], "Количество точек и векторов должно совпадать"


# Извлечение координат точек и компонентов векторов
X, Y = points[:, 0], points[:, 1]
U, V = vector_function(X, Y)

# Отрисовка векторного поля с уменьшенной толщиной линий и увеличенными векторами
plt.figure(figsize=(10, 10))
plt.quiver(X, Y, U, V, angles='xy', scale_units='xy', scale=0.9, width=0.004)
plt.xlim([X.min() - 1, X.max() + 1])
plt.ylim([Y.min() - 1, Y.max() + 1])
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Vector Field')
plt.grid()
plt.show()

