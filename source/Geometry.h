// Contributors Justin Hutchison yibbidy@gmail.com

#ifndef Geometry_h
#define Geometry_h

#include <math.h>
#include <list>
#include <vector>
#include <algorithm>

#define PI 3.141592654

// A vector class for doing vector math.
template<typename T, int Length> class TVector 
{
public:
	TVector() {
	}

	TVector(const TVector<T, Length> & vector) {
		*this = vector;
	}
	
	TVector(T x) {
		X() = x;			
	}
	
	TVector(T x, T y) {
		X() = x;
		Y() = y;
	}
	
	TVector(T x, T y, T z) {
		X() = x;
		Y() = y;
		Z() = z;
	}
	
	TVector(T x, T y, T z, T w) {
		X() = x;
		Y() = y;
		Z() = z;
		W() = w;
	}

	T & X() { return data[0]; }
	T & Y() { return data[1]; }
	T & Z() { return data[2]; }
	T & W() { return data[3]; }

	T & I() { return data[0]; }
	T & J() { return data[1]; }
	T & K() { return data[2]; }

	T & R() { return data[0]; }
	T & G() { return data[1]; }
	T & B() { return data[2]; }
	T & A() { return data[3]; }

	template<int SubLength> TVector<T, SubLength> ExtractSubVector(int begin) {
		TVector<T, SubLength> result;
		for( int i=begin; i<begin+SubLength; ++i ) result(i) = data[i];
		return result;
	}

	double Dot(TVector<T, Length> vector) {
		double result = 0.0;

		for( int i=0; i<Length; ++i ) result += data[i] * vector(i);

		return result;
	}

	TVector<T, Length> operator *(T uniformVector) {
		TVector<T, Length> result;

		for( int i=0; i<Length; ++i ) result(i) = data[i] * uniformVector;

		return result;
	}

	TVector<T, Length> operator /(T uniformVector) {
		TVector<T, Length> result;

		for( int i=0; i<Length; ++i ) result(i) = data[i] / uniformVector;

		return result;
	}

	TVector<T, Length> operator +(TVector<T, Length> vector) {
		TVector<T, Length> result;

		for( int i=0; i<Length; ++i ) result.data[i] = data[i] + vector(i);

		return result;
	}

	TVector<T, Length> operator +(T uniformVector) {
		TVector<T, Length> result;

		for( int i=0; i<Length; ++i ) result.data[i] = data[i]+uniformVector;

		return result;
	}

	TVector<T, Length> operator -(TVector<T, Length> vector) {
		TVector<T, Length> result;

		for( int i=0; i<Length; ++i ) result.data[i] = data[i]-vector(i);

		return result;
	}
	
	TVector<T, Length> operator -(T uniformVector) {
		TVector<T, Length> result;

		for( int i=0; i<Length; ++i ) result.data[i] = data[i]-uniformVector;

		return result;
	}
	
	TVector<T, Length> operator -() {
		TVector<T, Length> result;
		for( int i=0; i<Length; ++i ) result(i) = -data[i];
		return result;
	}
	
	TVector<T, Length> Cross(TVector<T, Length> matrix) {
		TVector<T, Length> result;
		result.X() = Y()*matrix.Z() - Z()*matrix.Y();
		result.Y() = -(X()*matrix.Z() - Z()*matrix.X());
		result.Z() = X()*matrix.Y() - Y()*matrix.X();
		return result;
	}

	double Magnitude() {
		double sum = 0.0;
		for( int i=0; i<Length; ++i ) sum += pow(data[i], 2);
		return sqrt(sum);
	}

	TVector<T, Length> Normalize() {
		TVector<T, Length> result;

		T magnitude = Magnitude();
		
		for( int i=0; i<Length; ++i )
			result(i) = data[i] / magnitude;
		
		return result;
	}
	
	double DistanceTo(TVector<T, Length> vector) {
		double sum = 0.0;
		for( int i=0; i<Length; ++i ) sum += pow(data[i]-vector(i), 2);
		return sqrt(sum);
	}

	double AngleBetween(TVector<T, Length> vector) {
		double cos = (Dot(vector)) / (Magnitude()*vector.Magnitude());
		if( cos<-1.0 || cos>1.0 ) return 0.0;
		else return acos(cos);
	}

	bool ParallelTo(TVector<T, Length> vector) {
		TVector<T, Length> v1 = vector.Normalize();
		TVector<T, Length> v2 = Normalize();
		
		if( (abs(X()-vector.X())<1e-4 && abs(Y()-vector.Y())<1e-4 && abs(Z()-vector.Z())<1e-4) ||
			(abs(X()+vector.X())<1e-4 && abs(Y()+vector.Y())<1e-4 && abs(Z()+vector.Z())<1e-4) ) return true;
		return false;
	}

	bool operator ==(TVector<T, Length> vector) {
		for( int i=0; i<Length; ++i ) {
			if( abs(vector(i) - data[i]) > 1e-5 ) return false;
		}
		return true;
	}
	
	bool operator !=(TVector<T, Length> vector) {
		return !(vector == *this);
	}
	
	T & operator ()(int i) {
		return data[i];
	}


	T data[Length];
};


// A matrix class for doing math.
template<typename T, int Rows, int Cols> class TMatrix 
{
public:
	operator TVector<T, Rows*Cols>() {
		TVector<T, Rows*Cols> vector;
		int i=0;
		for( int row=0; row<Rows; ++row ) {
			for( int col=0; col<Cols; ++col ) {
				vector(i++) = fData[row][col];
			}
		}
		return vector;
	}
	void operator =(const T data[Rows][Cols]) {
		for( int row=0; row<Rows; ++row ) {
			for( int col=0; col<Cols; ++col ) {
				this->fData[row][col] = data[row][col];
			}
		}
	}
	void operator =(const TMatrix<T, Rows, Cols> & matrix) {
		*this = matrix.fData;
	}

	static TMatrix<T, 1, Cols> FromRowVector(TVector<T, Cols> vector) {
		TMatrix<T, 1, Cols> matrix;
		for( int col=0; col<Cols; ++col ) matrix.fData[0][col] = vector(col);
		return matrix;
	}
	static TMatrix<T, Rows, 1> FromColVector(TVector<T, Rows> vector) {
		TMatrix<T, Rows, 1> matrix;
		for( int row=0; row<Rows; ++row ) matrix.fData[row][0] = vector(row);
		return matrix;
	}
	
	T & operator()(int row, int col) {
		return fData[row][col];
	}


	TMatrix() {
	}
	
	TMatrix(TMatrix<T, Rows, Cols> & matrix) {
		*this = matrix;
	}
	
	TMatrix(T data[Rows][Cols]) {
		*this = data;
	}
	
	TMatrix(T initialValue) {
		for( int row=0; row<Rows; ++row ) {
			for( int col=0; col<Cols; ++col ) {
				fData[row][col] = initialValue;
			}
		}
	}
		
	static TMatrix<T, Rows, Cols> FromEulerParameters(T e0, T e1, T e2, T e3) {
		TMatrix<T, Rows, Cols> matrix;
		matrix.LoadIdentity();

		matrix.fData[0][0] = pow(e0, 2)+pow(e1, 2)-pow(e2, 2)-pow(e3, 2);
		matrix.fData[1][0] = 2.0*(e1*e2+e0*e3);
		matrix.fData[2][0] = 2.0*(e1*e3-e0*e2);
		matrix.fData[0][1] = 2.0*(e1*e2-e0*e3);
		matrix.fData[1][1] = pow(e0, 2)+pow(e2, 2)-pow(e1, 2)-pow(e3, 2);
		matrix.fData[2][1] = 2.0*(e2*e3+e0*e1);
		matrix.fData[0][2] = 2.0*(e1*e3+e0*e2);
		matrix.fData[1][2] = 2.0*(e2*e3-e0*e1);
		matrix.fData[2][2] = pow(e0, 2)+pow(e3, 2)-pow(e1, 2)-pow(e2, 2);

		return matrix;
	}

	TMatrix<T, Rows, Cols> Transpose() {
		TMatrix<T, Rows, Cols> result;
		for( int row=0; row<Rows; ++row ) {
			for( int col=0; col<Cols; ++col ) {
				result(col, row) = fData[row][col];
			}
		}
		return result;
	}

	double Determinate3x3() {
		double determinate = 0.0;

		determinate =  fData[0][0] * (fData[1][1]*fData[2][2] - fData[1][2]*fData[2][1]);
		determinate += fData[0][1] * (fData[1][0]*fData[2][2] - fData[1][2]*fData[2][0]);
		determinate += fData[0][2] * (fData[1][0]*fData[2][1] - fData[1][1]*fData[2][0]);
		
		return determinate;
	}

	
	template<int OperandRows, int OperandCols>
	TMatrix<T, Rows, OperandCols> operator *(TMatrix<T, OperandRows, OperandCols> & matrix) {
		TMatrix<T, Rows, OperandCols> result;

		for( int row=0; row<Rows; ++row ) {
			for( int col=0; col<OperandCols; ++col ) {
				TMatrix<T, 1, Cols> rowMatrix = ExtractSubMatrix<1, Cols>(row, 0);
				TVector<T, Cols> rowVector = rowMatrix;
				TMatrix<T, OperandRows, 1> colMatrix = matrix.ExtractSubMatrix<OperandRows, 1>(0, col);
				TVector<T, OperandRows> colVector = colMatrix;

				result(row, col) = rowVector * colVector;
			}				
		}

		return result;
	}

	TVector<T, Rows> operator *(TVector<T, Rows> & vector) {
		TMatrix<T, Rows, 1> matrix = TMatrix<T, Rows, 1>::FromColVector(vector);
		return *this * matrix;
	}

	TMatrix<T, Rows, Cols> MultiplyElements(TMatrix<T, Rows, Cols> & matrix) {
		TMatrix<T, Rows, Cols> result;

		for( int row=0; row<Rows; ++row ) {
			for( int col=0; col<Cols; ++col ) {
				result(row, col) = fData[row][col] * matrix(row, col);
			}				
		}

		return result;
	}
	
	TMatrix<T, Rows, Cols> DivideElements(TMatrix<T, Rows, Cols> & matrix) {
		TMatrix<T, Rows, Cols> result;

		for( int row=0; row<Rows; ++row ) {
			for( int col=0; col<Cols; ++col ) {
				result(row, col) = fData[row][col] / matrix(row, col);
			}				
		}

		return result;
	}

	TMatrix<T, Rows, Cols> operator *(T d) {
		TMatrix<T, Rows, Cols> result;
		for( int row=0; row<Rows; ++row ) {
			for( int col=0; col<Cols; ++col ) {
				result(row, col) = fData[row][col] * d;
			}
		}
		return result;
	}
	
	TMatrix<T, Rows, Cols> operator /(T d) {
		TMatrix<T, Rows, Cols> result;
		for( int row=0; row<Rows; ++row ) {
			for( int col=0; col<Cols; ++col ) {
				result(row, col) = fData[row][col] / d;
			}
		}
		return result;
	}
	
	TMatrix<T, Rows, Cols> operator +(T d) {
		TMatrix<T, Rows, Cols> result;
		for( int row=0; row<Rows; ++row ) {
			for( int col=0; col<Cols; ++col ) {
				result(row, col) = fData[row][col] + d;
			}
		}
		return result;
	}
	
	TMatrix<T, Rows, Cols> operator -(T d) {
		TMatrix<T, Rows, Cols> result;
		for( int row=0; row<Rows; ++row ) {
			for( int col=0; col<Cols; ++col ) {
				result(row, col) = fData[row][col] - d;
			}
		}
		return result;
	}
		

	template<int SubMatrixRows, int SubMatrixCols>
	TMatrix<T, SubMatrixRows, SubMatrixCols> ExtractSubMatrix(int startRow, int startCol) {
		TMatrix<T, SubMatrixRows, SubMatrixCols> result;
		for( int row=startRow; row<startRow+SubMatrixRows; ++row ) {
			for( int col=startCol; col<startCol+SubMatrixCols; ++col ) {
				result(row-startRow, col-startCol) = fData[row][col];
			}
		}
		return result;
	}

	template<int SubMatrixRows, int SubMatrixCols>
	TMatrix<T, Rows, Cols> ReplaceSubMatrix(int startRow, int startCol, TMatrix<T, SubMatrixRows, SubMatrixCols> & matrix) {
		TMatrix<T, Rows, Cols> result = *this;

		for( int row=0; row<SubMatrixRows; ++row ) {
			for( int col=0; col<SubMatrixCols; ++col ) {
				result(row+startRow, col+startCol) = matrix(row, col);
			}
		}

		return result;
	}

	template<int SubMatrixRows, int SubMatrixCols>
	TMatrix<T, Rows+SubMatrixRows, Cols+SubMatrixCols> InsertSubMatrix(int startRow, int startCol, TMatrix<T, SubMatrixRows, SubMatrixCols> & matrix) {
		int rows = Rows+SubMatrixRows;
		int cols = Cols+SubMatrixCols;
		TMatrix<T, Rows+SubMatrixRows, Cols+SubMatrixCols> result;
		for( int row=0; row<rows; ++row ) {
			for( int col=0; col<cols; ++col ) {
				if( col<startCol ) {
					if( row<startRow ) result(row, col) = fData[row][col];
					else if( row<startRow+SubMatrixRows ) result(row, col) = 0;
					else result(row, col) = fData[row-SubMatrixRows][col];
				} else if( col<startCol+SubMatrixCols ) {
					if( row<startRow ) result(row, col) = 0;
					else if( row<startRow+SubMatrixRows ) result(row, col) = matrix(row-startRow, col-startCol);
					else result(row, col) = 0;
				} else {
					if( row<startRow ) result(row, col) = fData[row][col-SubMatrixCols];
					else if( row<startRow+SubMatrixRows ) result(row, col) = 0;
					else result(row, col) = fData[row-SubMatrixRows][col-SubMatrixCols];
				}
			}
		}
		return result;
	}

	
	TMatrix<T, Rows-1, Cols> RemoveRow(int rowToRemove) {
		TMatrix<T, Rows-1, Cols> result;
		int destRow=0;
		for( int row=0; row<Rows; ++row ) {
			if( row == rowToRemove ) continue;
			for( int col=0; col<Cols; ++col ) result(destRow, col) = fData[row][col];
			++destRow;
		}
		return result;
	}
	TMatrix<T, Rows, Cols-1> RemoveColumn(int columnToRemove) {
		TMatrix<T, Rows, Cols-1> result;
		int destColumn=0;
		for( int col=0; col<Cols; ++col ) {
			if( col == columnToRemove ) continue;
			for( int row=0; row<Rows; ++row ) result(row, destColumn) = fData[row][col];
			++destColumn;
		}
		return result;
	}
	
	static TMatrix<T, 4, 4> RotationMatrixAboutX(T angleInRadians) {
		TMatrix<T, 4, 4> rotation;
		rotation.LoadIdentity();
		rotation(0, 0) = cos(angleInRadians);
		rotation(0, 1) = sin(angleInRadians);
		rotation(1, 0) = -sin(angleInRadians);
		rotation(1, 1) = cos(angleInRadians);
		return rotation;
	}

	static TMatrix<T, 4, 4> RotationMatrixAboutY(T angleInRadians) {
		TMatrix<T, 4, 4> rotation;
		rotation.LoadIdentity();
		rotation(1, 1) = cos(angleInRadians);
		rotation(1, 2) = sin(angleInRadians);
		rotation(2, 1) = -sin(angleInRadians);
		rotation(2, 2) = cos(angleInRadians);
		return rotation;
	}

	static TMatrix<T, 4, 4> RotationMatrixAboutZ(T angleInRadians) {
		TMatrix<T, 4, 4> rotation;
		rotation.LoadIdentity();
		rotation(0, 0) = cos(angleInRadians);
		rotation(0, 2) = -sin(angleInRadians);
		rotation(2, 0) = sin(angleInRadians);
		rotation(2, 2) = cos(angleInRadians);
		return rotation;
	}

	static TMatrix<T, 4, 4> TranslationMatrix(TVector<T, 3> & translateVector) {
		TMatrix<T, 4, 4> translation;
		translation.LoadIdentity();
		translation(3, 0) = translateVector.X();
		translation(3, 1) = translateVector.Y();
		translation(3, 2) = translateVector.Z();
		return translation;
	}

	static TMatrix<T, 4, 4> ScaleMatrix(T scale) {
		TMatrix<T, 4, 4> result;
		result.LoadIdentity();
		result(0, 0) = result(1, 1) = result(2, 2) = scale;
		return result;			 
	}

	static TMatrix<T, 3, 3> ScaleMatrix(T xScale, T yScale) {
		TMatrix<T, 3, 3> result;
		result.LoadIdentity();
		result(0, 0) = xScale;
		result(1, 1) = yScale;
		return result;			 
	}

	TMatrix<T, Rows, Cols> RotationMatrixBetween(TMatrix<T, Rows, Cols> & matrix) {
		return Transpose() * matrix;
	}

	bool operator ==(TMatrix<T, Rows, Cols> & matrix) {
		for( int row=0; row<Rows; ++row ) {
			for( int col=0; col<Cols; ++col ) {
				if( abs(fData[row][col] - matrix(row, col)) > 1e-5 ) return false;
			}
		}

		return true;
	}
	bool operator !=(TMatrix<T, Rows, Cols> & matrix) {
		return !(*this == matrix);
	}

	void LoadZero() {
		for( int row=0; row<Rows; ++row ) {
			for( int col=0; col<Cols; ++col ) {
				fData[row][col] = 0;
			}
		}
	}

	void LoadIdentity() {
		if( Rows != Cols ) return;

		for( int row=0; row<Rows; ++row ) {
			for( int col=0; col<Cols; ++col ) {
				fData[row][col] = (row==col) ? 1 : 0;
			}
		}
	}


	TVector<T, Rows> Solve(TVector<T, Rows> rhs) {
		bool solvable;
		return Solve(rhs, solvable);
	}

	TVector<T, Rows> Solve(TVector<T, Rows> rhs, bool & solvable) {

		solvable = true;
		TMatrix<T, Rows, Cols> lhs = *this;

		// initialize reindex, solution
		TVector<T, Rows> solution;
		bool fixed[Rows];
		int reindex[Rows];
		for( int i=0; i<Rows; ++i ) {
			reindex[i] = i;
			solution(i) = 1;
			fixed[i] = false;
		}

		// organize the lhs
		for( int col=0; col<Cols; ++col ) {
			for( int row=0; row<Rows; ++row ) {
				if( !fixed[row] && lhs(row, col)!=0.0 ) {
					if( reindex[col] != row ) std::swap(reindex[col], reindex[row]);
					fixed[row] = true;

					T coeff = lhs(row, col);

					for( int c=0; c<Cols; ++c ) lhs(row, c) /= coeff;
					rhs(row) /= coeff;

					for( int r=0; r<Rows; ++r ) {
						if( r==row || 0.0==lhs(r, col) ) continue;
						coeff = lhs(r, col);
						for( int c=0; c<Cols; ++c) {
							lhs(r, c) -= coeff * lhs(row, c);
						}
						rhs(r) -= coeff * rhs(row);
					}
					break;
				}
			}
		}

		// case 1: zero column
		for( int i=0; i<Rows; ++i ) {
			if( !fixed[i] ) solution(reindex[i]) = 0.0;
		}

		// solve
		for( int row=0; row<Rows; ++row ) {
			if( 1.0 == lhs(reindex[row], row) ) {
				bool zeros = true;
				for( int col=0; col<Cols && zeros; ++col ) {
					if( col == row ) continue;
					if( lhs(reindex[row], col) != 0.0 ) zeros = false;
				}
				if( zeros ) {
					// case 3: easy case
					solution(reindex[row]) = rhs(row);
				} else {
					// case 4: hard case
					T sum = rhs(row);
					for( int col=0; col<Cols; ++col ) {
						if( col == row ) continue;
						sum -= lhs(reindex[row], col) * solution(col);
					}
					solution(reindex[row]) = sum;
				}
			} else {
				// case 2: should be all zeros
				for( int col=0; col<Cols; ++col ) {
					if( lhs(reindex[row], col) != 0.0 ) {
						solvable = false;
						break;
					}
					if( rhs(reindex[row]) != 0.0 ) solvable = false;
					if( !solvable ) return /*input error*/ solution;
				}
			}
		}
			
		return solution;
	}

	T fData[Rows][Cols];
};

template<typename T, int Rows, int Cols> TVector<T, Rows> operator *(TVector<T, Rows> vector, TMatrix<T, Rows, Cols> matrix) 
{
	return (TVector<T, Rows>) (TMatrix<T, 1, Rows>::FromRowVector(vector) * matrix);
}

// geometry with lines
class TLine 
{
public:

	TLine(TVector<double, 2> p0, TVector<double, 2> direction, bool segment) {
		this->fP0 = p0;
		this->fDirection = direction;
		this->fSegment = segment;
	}

	TVector<double, 2> Intersect(TLine & line, bool & intersect) {
		if( fDirection.ParallelTo(line.fDirection) ) {
			intersect = true;
			return fP0;
		}

		TMatrix<double, 2, 2> lhs;
		lhs(0, 0) = fDirection.X(); lhs(0, 1) = -line.fDirection.X();
		lhs(1, 0) = fDirection.Y(); lhs(1, 1) = -line.fDirection.Y();

		TVector<double, 2> rhs = line.fP0 - fP0;

		TVector<double, 2> solution = lhs.Solve(rhs, intersect);

		if( !intersect ) return fP0;

		TVector<double, 2> intersection = fP0 + fDirection*solution(0);

		if( fSegment && intersect ) {
            double distance0 = intersection.DistanceTo(fP0);
			double distance1 = intersection.DistanceTo(fP0+fDirection);
			double segmentLength = fDirection.Magnitude();
			if( distance0>segmentLength || distance1>segmentLength ) intersect = false;
		}

		if( line.fSegment && intersect ) {
            double distance0 = intersection.DistanceTo(line.fP0);
			double distance1 = intersection.DistanceTo(line.fP0+line.fDirection);
			double segmentLength = line.fDirection.Magnitude();
			if( distance0>segmentLength || distance1>segmentLength ) intersect = false;
		}

        return intersection;
	}

	TVector<double, 2> fP0;
	TVector<double, 2> fDirection;
	bool fSegment;
};

// geometry with polygons
class TPolygon 
{
public:
	TPolygon() {
	}

	TPolygon(std::list<TVector<double, 2> > vertices) {
		this->fVertices = vertices;
	}

	TVector<double, 2> FirstVertex(std::vector<TVector<double, 2> > & vertices1, std::vector<TVector<double, 2> > & vertices2, int & firstVertex, bool & disjoint) {
		disjoint = true;
		firstVertex = 0;

		if( vertices1.size()<=3 || vertices2.size()<=3 ) {
			disjoint = false;
			return TVector<double, 2>(0, 0);
		}

		bool vertexIsInPolygon = false;
		for( std::vector<TVector<double, 2> >::iterator it = vertices2.begin(); it != vertices2.end(); ++it ) {
			if( PointIsInPolygon(*it) ) {
				disjoint = false;
				return *it;
			}
			++firstVertex;
		}

		for( int i=0; i<(int)vertices1.size()-1; ++i ) {
			TVector<double, 2> p1 = vertices1[i];
			TVector<double, 2> p2 = vertices1[i+1];
			TLine l1(p1, p2-p1, true);
			for( int j=0; j<(int)vertices2.size()-1; ++j ) {
				TVector<double, 2> p3 = vertices2[i];
				TVector<double, 2> p4 = vertices2[i+1];
				TLine l2(p3, p4-p3, true);
				bool intersect;
				TVector<double, 2> firstVertex = l1.Intersect(l2, intersect);
				if( intersect ) {
					disjoint = false;
					firstVertex = i;
					return firstVertex; 
				}
			}
		}

		firstVertex = 0;
		return TVector<double, 2>(0, 0);
	}

	TPolygon Clip(TPolygon polygon) {
		std::vector<TVector<double, 2> > vertices1(fVertices.size()); copy(fVertices.begin(), fVertices.end(), vertices1.begin());
		std::vector<TVector<double, 2> > vertices2(polygon.fVertices.size()); copy(polygon.fVertices.begin(), polygon.fVertices.end(), vertices2.begin());

		std::vector<TVector<double, 2> > * polyVertices[] = { &vertices1, &vertices2 };
		int vertexCount[] = { (int)fVertices.size(), (int)polygon.fVertices.size() };
		
		int currentPoly = 1;
		int currentVertex = 0;
		int lastEdge = -1;

		bool disjoint;
		TVector<double, 2> firstVertex = FirstVertex(vertices1, vertices2, currentVertex, disjoint);
		if( disjoint ) return TPolygon();

		std::list<TVector<double, 2> > clipVertices;
		clipVertices.push_back(firstVertex);

		int loopCount = 0;
		do {
			int nextPoly = (currentPoly+1) % 2;

			TVector<double, 2> p1 = (*polyVertices[currentPoly])[currentVertex];
			TVector<double, 2> p2 = (*polyVertices[currentPoly])[(currentVertex+1) % vertexCount[currentPoly]];
			TLine l1(p1, p2-p1, true);

			bool intersect = false;

			for( int edge=0; edge<vertexCount[nextPoly]; ++edge ) {
				if( edge == lastEdge ) continue;

				TVector<double, 2> p3 = (*polyVertices[nextPoly])[edge];
				TVector<double, 2> p4 = (*polyVertices[nextPoly])[(edge+1) % vertexCount[nextPoly]];
				TLine l2(p3, p4-p3, true);

				TVector<double, 2> intersection = l1.Intersect(l2, intersect);
				if( intersect ) {
					if( find(clipVertices.begin(), clipVertices.end(), intersection) != clipVertices.end() ) {
						// TODO is this correct? no it isn't
						return TPolygon(clipVertices);
					}
					clipVertices.push_back(intersection);
					lastEdge = currentVertex;
					currentPoly = (currentPoly+1) % 2;
					currentVertex = edge;
					break;
				}
				
			}

			if( !intersect ) {
				currentVertex = (currentVertex+1) % vertexCount[currentPoly];
				clipVertices.push_back((*polyVertices[currentPoly])[currentVertex]);
				lastEdge = -1;
			}

			if( ++loopCount >= 100 ) break;
		} while( clipVertices.back() != firstVertex );
		
		clipVertices.pop_back();

		return TPolygon(clipVertices);
	}

	bool PointIsInPolygon(TVector<double, 2> point) {
		double angleSum = 0.0;

		TVector<double, 2> previous = fVertices.back();

		for( std::list<TVector<double, 2> >::iterator it = fVertices.begin(); it != fVertices.end(); ++it ) {
			TVector<double, 2> current = *it;

			TVector<double, 2> v1 = previous-point;
			TVector<double, 2> v2 = current-point;

			// perpendicular vector is clockwise from v1
			TVector<double, 2> perpendicular(v1.Y(), -v1.X());
			double sign = ((perpendicular.Dot(v2)) > 0) ? 1.0 : -1.0;

			angleSum += sign*v1.AngleBetween(v2);

			previous = current;
		}

		double epsilon = 1e-3;
		return abs(abs(angleSum) - 2*PI) < epsilon;
	}

	std::list<TVector<double, 2> > fVertices;

};

// geometry with rectangles
class TRectangle 
{
public:
	TRectangle() {
	};

	TRectangle(TVector<double, 2> & bottomLeft, double width, double height) {
		this->fBottomLeft = bottomLeft;
		this->fWidth = width;
		this->fHeight = height;
	}

	bool ContainsPoint(TVector<double, 2> & point) {
		double left = fBottomLeft.X();
		double right = left + fWidth;
		double bottom = fBottomLeft.Y();
		double top = bottom + fHeight;

		return !(point.X()<left || point.X()>right || point.Y()<bottom || point.Y()>top);
	}

	bool Intersect(TRectangle & rectangle) {
		double left1 = fBottomLeft.X();
		double right1 = left1 + fWidth;
		double bottom1 = fBottomLeft.Y();
		double top1 = bottom1 + fHeight;

		double left2 = rectangle.fBottomLeft.X();
		double right2 = left2 + rectangle.fWidth;
		double bottom2 = rectangle.fBottomLeft.Y();
		double top2 = bottom2 + rectangle.fHeight;

		return !(left1>right2 || right1<left2 || bottom1>top2 || top1<bottom2);
	}

	TVector<double, 2> Center() {
		return fBottomLeft + TVector<double, 2>(fWidth, fHeight) / 2;
	}

	TVector<double, 2> fBottomLeft;
	double fWidth, fHeight;
};

#endif

