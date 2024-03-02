#include <algorithm>
#include <cstdint>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>

// Класс матриц
template <size_t N, size_t M, typename T = int64_t>
class Matrix {
 public:
  Matrix() { matr_.resize(N, std::vector<T>(M)); }

  Matrix(const std::vector<std::vector<T>>& from) { matr_ = from; }

  Matrix(const T& elem) {
    for (size_t i = 0; i < N; ++i) {
      std::vector<T> temp;
      for (size_t j = 0; j < M; ++j) {
        temp.push_back(elem);
      }
      matr_.push_back(temp);
    }
  }

  Matrix(const Matrix& from) = default;

  Matrix<N, M, T>& operator=(const Matrix& from) = default;

  ~Matrix() = default;

  Matrix<N, M, T>& operator+=(const Matrix<N, M, T>& add) {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        matr_[i][j] += add(i, j);
      }
    }
    return *this;
  }

  Matrix<N, M, T>& operator-=(const Matrix<N, M, T>& dec) {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        matr_[i][j] -= dec(i, j);
      }
    }
    return *this;
  }

  Matrix<N, M, T>& operator*=(const T& elem) {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        matr_[i][j] *= elem;
      }
    }
    return *this;
  }

  template <size_t K>
  Matrix<N, K, T> operator*(Matrix<M, K, T> mnoz) {
    Matrix<N, K, T> ans;
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < K; ++j) {
        for (size_t iter = 0; iter < M; ++iter) {
          ans(i, j) += matr_[i][iter] * mnoz(iter, j);
        }
      }
    }
    return ans;
  }

  Matrix<N, M, T> operator+(const Matrix<N, M, T>& add1) {
    return (Matrix<N, M, T>(*this) += add1);
  }

  Matrix<N, M, T> operator-(const Matrix<N, M, T>& dec1) {
    return (Matrix<N, M, T>(*this) -= dec1);
  }

  Matrix<N, M, T> operator*(const T& elem) {
    return (Matrix<N, M, T>(*this) *= elem);
  }

  Matrix<M, N, T> Transposed() {
    Matrix<M, N, T> temp;
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        temp(j, i) = matr_[i][j];
      }
    }
    return temp;
  }

  T& operator()(const size_t& str, const size_t& stl) {
    return matr_[str][stl];
  }

  const T& operator()(const size_t& str, const size_t& stl) const {
    return matr_[str][stl];
  }

  bool operator==(const Matrix<N, M, T>& mat) const {
    return (matr_ == mat.matr_);
  }

  bool operator!=(const Matrix<N, M, T>& mat) const {
    return !(matr_ == mat.matr_);
  }

  void SwapRows(const size_t& from, const size_t& to) {
    std::swap(matr_[from], matr_[to]);
  }

  void PrintMatrix() {
    for (size_t indx = 0; indx < N; ++indx) {
      for (size_t jndx = 0; jndx < M; ++jndx) {
        std::cout << matr_[indx][jndx] << ' ';
      }
      std::cout << '\n';
    }
  }

 private:
  std::vector<std::vector<T>> matr_;
};

// Функция получения данных из файла
template <size_t N, size_t M, typename T = int64_t>
void GetData(std::string file_name, Matrix<N, M, T>& first_matrix,
             Matrix<N, 1, T>& second_matrix) {
  // окрываем файл для чтения
  std::ifstream in(file_name);
  if (in.is_open()) {
    size_t indx = 0;
    size_t jndx = 0;
    float nmb_in;
    std::vector<std::vector<float>> in_vector(N);
    for (size_t kndx = 0; kndx < N; ++kndx) {
      in_vector[kndx].resize(M);
    }
    // считываем данные в большую матрицу
    while ((indx != N) && (jndx != M)) {
      in >> nmb_in;
      in_vector[indx][jndx] = nmb_in;
      ++jndx;
      if (jndx == M) {
        jndx = 0;
        ++indx;
      }
    }

    first_matrix = in_vector;
    indx = 0;
    jndx = 0;
    for (size_t kndx = 0; kndx < N; ++kndx) {
      in_vector[kndx].resize(1);
    }
    // считываем данные в маленькую матрицу
    while (indx != N) {
      in >> nmb_in;
      in_vector[indx][jndx] = nmb_in;
      ++indx;
    }
    second_matrix = in_vector;
  }
  // закрываем файл для чтения
  in.close();
}

// LUP декомпозиция, возвращает матрицу С вида C = L + U - E
template <size_t N, size_t M, typename T = int64_t>
Matrix<N, M, T> LUPDecompose(Matrix<N, M, T>& A, Matrix<N, M, T>& P) {
  Matrix<N, M, T> C = A;
  for(size_t indx = 0; indx < N; ++indx) {
    //поиск опорного элемента
    float pivotValue = -1;
    size_t pivot = SIZE_MAX;
    for(size_t row = indx; row < N; ++row) {
      if(std::fabs(C(row, indx)) > pivotValue) {
        pivotValue = std::fabs(C(row, indx));
        pivot = row;
      }
    }

    if(pivotValue != -1) {
      //меняем местами i-ю строку и строку с опорным элементом
      P.SwapRows(pivot, indx);
      C.SwapRows(pivot, indx);

      // делим элементы ниже i на сам опорный элемент
      for(size_t jndx = indx + 1; jndx < N; ++jndx) {
        C(jndx, indx) /= C(indx, indx);
        for(size_t kndx = indx + 1; kndx < N; ++kndx)
          C(jndx, kndx) -= C(jndx , indx) * C(indx, kndx);
      }

    }
  }
  return C;
}

// Основная функция решения СЛАУ
template <size_t N, size_t M, typename T = int64_t>
Matrix<N, 1, T> Solve(Matrix<N, M, T>& A, Matrix<N, 1, T>& f) {
  std::vector<std::vector<float>> identity_matrix(N);
  for (size_t indx = 0; indx < N; ++indx) {
    identity_matrix[indx].resize(M);
    identity_matrix[indx][indx] = 1;
  }

  // обозначаем Р как единичную матрицу
  Matrix<N, M, T> P = identity_matrix;

  // находим разложение матрицы А
  Matrix<N, M, T> C = LUPDecompose(A, P);

  // преобразуем правую часть
  Matrix<N, 1, float> f_new = P * f;

  // вычисляем верхнюю и нижнюю треугольные матрицы используя свойства матрицы С
  Matrix<N, M, T> L = C;
  for (size_t indx = 0; indx < N; ++indx) {
    for (size_t jndx = indx; jndx < M; ++jndx) {
      L(indx, jndx) = 0;
    }
  }
  Matrix<N, M, T> U = C - L;
  for (size_t indx = 0; indx < N; ++indx) {
    L(indx, indx) = 1;
  }

  // собственно решаем уравнение
  std::vector<float> y(N);
  float sum;
  int64_t indx;
  int64_t jndx;

  // делаем прямой ход и обратный ход
  for(indx = 0; indx < N; ++indx) {
    sum = 0;
    for(jndx = 0; jndx <= indx - 1; ++jndx)
      sum += L(indx, jndx) * y[jndx];
    y[jndx] = (f_new(indx, 0) - sum) / L(indx, indx);
  }

  // находим вектор решения
  Matrix<N, 1, float> x(0);
  for(indx = N - 1; indx >=0 ; --indx) {
    sum = 0;
    for(jndx = N - 1; jndx > indx; --jndx)
      sum += U(indx, jndx) * x(jndx, 0);
    x(indx, 0) = (y[indx] - sum) / U(indx, indx);
  }

  std::cout << '\n' << "Matrix P:" << '\n';
  P.PrintMatrix();
  std::cout << '\n' << "Matrix U:" << '\n';
  U.PrintMatrix();
  std::cout << '\n' << "Matrix L:" << '\n';
  L.PrintMatrix();
  std::cout << '\n' << "Matrix for answer:" << '\n';
  x.PrintMatrix();
  return x;
}

// Простейшая функция нормы двух векторов
template <size_t N, size_t M, typename T = int64_t>
float Norma(Matrix<N, M, T>& first, Matrix<N, M, T>& second) {
  Matrix<N, M, T> matrix = first - second;
  for (size_t indx = 0; indx < N; ++indx) {
    matrix(indx, 0) = std::pow(matrix(indx, 0), 2);
  }
  float sum = 0;
  for (size_t indx = 0; indx < N; ++indx) {
    sum += matrix(indx, 0);
  }
  return sqrtf(sum);
}

int main() {
  Matrix<3, 3, float> A1;
  Matrix<4, 4, float> A2;
  Matrix<3, 3, float> A3;
  Matrix<10, 10, float> A4;
  Matrix<3, 1, float> B1;
  Matrix<4, 1, float> B2;
  Matrix<3, 1, float> B3;
  Matrix<10, 1, float> B4;

  // получаем данные матриц из файлов с ними
  GetData("../A1 B1", A1, B1);
  GetData("../A2 B2", A2, B2);
  GetData("../A3 B3", A3, B3);
  GetData("../A4 B4", A4, B4);

  // находим решения с помощью LUP разложения и находим
  // разницу между эталонным решением и получившимся

  // решаем первый пример
  std::cout << '\n' << "Solutions for first matrices:" << '\n';
  auto x1 = Solve(A1, B1);
  std::vector<std::vector<float>> tmp = {{-8.17310}, {5.15291}, {3.47127}};
  Matrix<3, 1, float> known_ans_1 = tmp;
  std::cout << '\n' << "The difference between the received solution and the reference one: "
            << Norma(x1, known_ans_1) << '\n';

  // решаем второй пример
  std::cout << '\n' << "Solutions for second matrices:" << '\n';
  auto x2 = Solve(A2, B2);
  tmp = {{0.40970}, {-0.32679}, {-0.45735}, {0.54471}};
  Matrix<4, 1, float> known_ans_2 = tmp;
  std::cout << '\n' << "The difference between the received solution and the reference one: "
            << Norma(x2, known_ans_2) << '\n';

  // решаем третий пример
  std::cout << '\n' << "Solutions for third matrices:" << '\n';
  auto x3 = Solve(A3, B3);
  tmp = {{0.54277}, {1.04385} , {1.57772}};
  Matrix<3, 1, float> known_ans_3 = tmp;
  std::cout << '\n' << "The difference between the received solution and the reference one: "
            << Norma(x3, known_ans_3) << '\n';

  // решаем четвертый пример
  std::cout << '\n' << "Solutions for fourth matrices:" << '\n';
  auto x4 = Solve(A4, B4);
  tmp = {{0.70182}, {-0.13286} , {0.39599}, {0.22257}, {0.27678},
         {-0.06892}, {-0.66132}, {-1.39346}, {-0.40983}, {-0.74357}};
  Matrix<10, 1, float> known_ans_4 = tmp;
  std::cout << '\n' << "The difference between the received solution and the reference one: "
            << Norma(x4, known_ans_4) << '\n';

  return 0;
}
