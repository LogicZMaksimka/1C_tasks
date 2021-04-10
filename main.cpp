#include <iostream>
#include <vector>

using MatrixType = int;

template<typename T>
T Abs(const T &x) {
  return (x < 0) ? -x : x;
}

std::vector<std::vector<MatrixType>> DeleteCross(const std::vector<std::vector<
	MatrixType>> &matrix,
												 int i, int j) {
  int matrix_size = matrix.size();
  std::vector<std::vector<MatrixType>>
	  new_matrix(matrix_size - 1, std::vector<MatrixType>(matrix_size - 1));
  for (int k = 0; k < matrix_size - 1; ++k) {
	for (int l = 0; l < matrix_size - 1; ++l) {
	  new_matrix[k][l] = matrix[k + ((k >= i) ? 1 : 0)][l + ((l >= j) ? 1 : 0)];
	}
  }
  return new_matrix;
}

int GetDeterminant(const std::vector<std::vector<MatrixType>> &matrix) {
  if (matrix.size() == 1) {
	return matrix[0][0];
  }
  if (matrix.size() == 2) {
	return matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];
  }
  int det = 0;
  for (int i = 0; i < matrix.size(); ++i) {
	det += ((i % 2 == 0) ? 1 : -1) * matrix[0][i]
		* GetDeterminant(DeleteCross(matrix, 0, i));
  }
  return det;
}

enum DiagonalType {
  kMainDiagonal,
  kSideDiagonal
};

std::vector<std::vector<MatrixType>> TransposeMatrix(
	const std::vector<std::vector<MatrixType>> &matrix,
	DiagonalType diagonal_type) {
  int matrix_size = matrix.size();
  std::vector<std::vector<MatrixType>>
	  new_matrix(matrix_size, std::vector<MatrixType>(matrix_size));
  if (diagonal_type == kMainDiagonal) {
	for (int i = 0; i < matrix_size; ++i) {
	  for (int j = 0; j < matrix_size; ++j) {
		new_matrix[j][i] = matrix[i][j];
	  }
	}
  } else {
	for (int i = 0; i < matrix_size; ++i) {
	  for (int j = 0; j < matrix_size; ++j) {
		new_matrix[i][j] = matrix[matrix_size - j - 1][matrix_size - i - 1];
	  }
	}
  }
  return new_matrix;
}

std::vector<std::vector<MatrixType>> TransposeSubmatrix(
	std::vector<std::vector<MatrixType>> matrix,
	DiagonalType diagonal_type,
	int upper_left_x,
	int upper_left_y,
	int submatrix_size) {

  std::vector<std::vector<MatrixType>>
	  submatrix(submatrix_size, std::vector<MatrixType>(submatrix_size));

  int lower_right_x = upper_left_x + submatrix_size;
  int lower_right_y = upper_left_y + submatrix_size;
  for (int i = upper_left_y; i < lower_right_y; ++i) {
	for (int j = upper_left_x; j < lower_right_x; ++j) {
	  submatrix[i - upper_left_y][j - upper_left_x] = matrix[i][j];
	}
  }

  std::vector<std::vector<MatrixType>>
	  transposed_submatrix = TransposeMatrix(submatrix, diagonal_type);

  for (int i = upper_left_y; i < lower_right_y; ++i) {
	for (int j = upper_left_x; j < lower_right_x; ++j) {
	  matrix[i][j] = transposed_submatrix[i - upper_left_y][j - upper_left_x];
	}
  }

  return matrix;
}

// полное переборное решение:
// транспонирует каждую подматрицу в матрице
// и входит в рекурсию относительно новой матрицы

// этот перебор не самый эффективный: лучше было бы делать
// преобразования матриц по уровням - так чтобы на каждом уровне
// были все матрицы с предыдущего уровня, к которым применили всевозможные
// преобразования подматриц. После этого можно было удалить совпадающие матрицы
// и сильно уменьшить время работы
std::vector<std::vector<MatrixType>> GetMatrixWithMaxDeterminant(
	const std::vector<std::vector<MatrixType>> &matrix, int depth) {
  if (depth == 0) {
	return matrix;
  }

  int matrix_size = matrix.size();
  int max_det = Abs(GetDeterminant(matrix));
  std::vector<std::vector<MatrixType>> matrix_with_max_det = matrix;
  for (int submatrix_size = 2; submatrix_size < matrix_size; ++submatrix_size) {
	for (int submatrix_x_coord = 0;
		 submatrix_x_coord < matrix_size - submatrix_size + 1;
		 ++submatrix_x_coord) {
	  for (int submatrix_y_coord = 0;
		   submatrix_y_coord < matrix_size - submatrix_size + 1;
		   ++submatrix_y_coord) {

		std::vector<std::vector<MatrixType>>
			new_submatrix_main_diag = TransposeSubmatrix(matrix,
														 kMainDiagonal,
														 submatrix_x_coord,
														 submatrix_y_coord,
														 submatrix_size);

		std::vector<std::vector<MatrixType>>
			new_submatrix_side_diag = TransposeSubmatrix(matrix,
														 kSideDiagonal,
														 submatrix_x_coord,
														 submatrix_y_coord,
														 submatrix_size);
		int main_det = Abs(GetDeterminant(GetMatrixWithMaxDeterminant(
			new_submatrix_main_diag, depth - 1)));
		int side_det = Abs(GetDeterminant(GetMatrixWithMaxDeterminant(
			new_submatrix_side_diag, depth - 1)));

		if (main_det > max_det) {
		  matrix_with_max_det = new_submatrix_main_diag;
		  max_det = main_det;
		}
		if (side_det > max_det) {
		  matrix_with_max_det = new_submatrix_side_diag;
		  max_det = side_det;
		}
	  }
	}
  }
  return matrix_with_max_det;
}

// формат ввода:
// размер матрицы
// сама матрица

int main() {
  int matrix_size;
  std::cin >> matrix_size;
  std::vector<std::vector<MatrixType>>
	  matrix(matrix_size, std::vector<MatrixType>(matrix_size));

  for (int i = 0; i < matrix_size; ++i) {
	for (int j = 0; j < matrix_size; ++j) {
	  std::cin >> matrix[i][j];
	}
  }

  // можно рассчитать точную формулу глубины рекурсии
  // на которой будут перебраны все возможные комбинации
  int depth = 1;
  std::vector<std::vector<MatrixType>>
	  new_matrix = GetMatrixWithMaxDeterminant(matrix, depth);


  std::cout << "Искомая матрица: \n";
  for (int i = 0; i < new_matrix.size(); ++i) {
	for (int j = 0; j < new_matrix.size(); ++j) {
	  std::cout << new_matrix[i][j] << " ";
	}
	std::cout << std::endl;
  }
  std::cout << "детерминант: " << GetDeterminant(new_matrix);
  return 0;
}


