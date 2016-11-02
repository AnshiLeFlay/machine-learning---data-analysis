function kernels(u, v, type) {
	if (type == 0) { //<u, v>^2, X=R^2
		var ans = Math.pow(u[0], 2) * Math.pow(v[0], 2) + Math.pow(u[1], 2) * Math.pow(v[1], 2) + 2 * u[0] * v[0] * u[1] * v[1];
	}
	return ans;
}

//Cholesky decomposition
function cholesky(Qss) {
	var n = Qss.length;
	var ret = [];

	for (var r = 0; r < n; r++) {
		ret[r] = [];
	    for (var c = 0; c <= r; c++)
	    {
	        if (c == r)
	        {
	            var sum = 0;
	            for (var j = 0; j < c; j++)
	            {
	                sum += ret[c][j] * ret[c][j];
	            }
	            ret[c][c] = Math.sqrt(Qss[c][c] - sum);
	        }
	        else
	        {
	            var sum = 0;
	            for (var j = 0; j < c; j++)
	                sum += ret[r][j] * ret[c][j];
	            ret[r][c] = 1 / ret[c][c] * (Qss[r][c] - sum);
	        }
	    }
	}

	return ret;
}

function transp(mtrx) {
	var ans = [];
	var m = mtrx.length;
	var n = mtrx[0].length;
	for (var i = 0; i < n; i++) {
		ans[i] = [];
		for (var j = 0; j < m; j++) {
			ans[i][j] = mtrx[j][i];
		}
	}
	return ans;
}

function mult(mtrxA, mtrxB) {
	var ans = [];

	var l = mtrxA.length;
	var m = mtrxA[0].length;

	var mB = mtrxB.length;
	var n = mtrxB[0].length;

	if (m == mb) {
		for (var i = 0; i < l; i++) {
			ans[i] = [];
			for (var j = 0; j < n; j++) {
				ans[i][j] = 0;
				for (var r = 0; r < m; r++) {
					ans[i][j] += mtrxA[i][r] * mtrxB[r][j];
				}
			}
		}
		return ans;
	}
}

/*
function detG(mtrx) {
	var ans; 

}

function detT(mtrx) {
	var ans = 1;

	for (var i = 0; i < mtrx.length; i++) {
		ans = ans * mtrx[i][i];
	}

	return ans;
}
function minorT(mtrx, n, m) {
	var ans = 1;

	for (var i = 0; i < mtrx.length; i++) {
		if (i != n)
			ans = ans * mtrx[i][i];
	}

	return ans;
}

function invert(mtrx) {
	var ans = [];
	var det = detT(mtrx);
	//var minor;

	for (var i = 0; i < mtrx.length; i++) {
		ans[i] = [];
		for (var j = 0; j < mtrx[0].length) {
			ans[i][j] = Math.pow(-1, i+j) * minorT(mtrx, i, j);
		}
	}

	return ans;
}
*/
//////////////////////////////////////////////////////////////////////////////////////////////
//will be changed with Gauss method
//////////////////////////////////////////////////////////////////////////////////////////////
function Determinant(A)
{
    var N = A.length, B = [], denom = 1, exchanges = 0;
    for (var i = 0; i < N; ++i)
     { B[i] = [];
       for (var j = 0; j < N; ++j) B[i][j] = A[i][j];
     }
    for (var i = 0; i < N-1; ++i)
     { var maxN = i, maxValue = Math.abs(B[i][i]);
       for (var j = i+1; j < N; ++j)
        { var value = Math.abs(B[j][i]);
          if (value > maxValue){ maxN = j; maxValue = value; }
        }
       if (maxN > i)
        { var temp = B[i]; B[i] = B[maxN]; B[maxN] = temp;
          ++exchanges;
        }
       else { if (maxValue == 0) return maxValue; }
       var value1 = B[i][i];
       for (var j = i+1; j < N; ++j)
        { var value2 = B[j][i];
          B[j][i] = 0;
          for (var k = i+1; k < N; ++k) B[j][k] = (B[j][k]*value1-B[i][k]*value2)/denom;
        }
       denom = value1;
     }
    if (exchanges%2) return -B[N-1][N-1];
    else return B[N-1][N-1];
}
function AdjugateMatrix(A)
{                                        
    var N = A.length, adjA = [];
    for (var i = 0; i < N; i++)
     { adjA[i] = [];
       for (var j = 0; j < N; j++)
        { var B = [], sign = ((i+j)%2==0) ? 1 : -1;
          for (var m = 0; m < j; m++)
           { B[m] = [];
             for (var n = 0; n < i; n++)   B[m][n] = A[m][n];
             for (var n = i+1; n < N; n++) B[m][n-1] = A[m][n];
           }
          for (var m = j+1; m < N; m++)
           { B[m-1] = [];
             for (var n = 0; n < i; n++)   B[m-1][n] = A[m][n];
             for (var n = i+1; n < N; n++) B[m-1][n-1] = A[m][n];
           }
          adjA[i][j] = sign*Determinant(B);
        }
     }
    return adjA;
}
function InverseMatrix(A)
{   
    var det = Determinant(A);
    if (det == 0) return false;
    var N = A.length, A = AdjugateMatrix(A);
    for (var i = 0; i < N; i++)
     { for (var j = 0; j < N; j++) A[i][j] /= det; }
    return A;
}
//////////////////////////////////////////////////////////////////////////////////////////////

var Ls = cholesky(Qss);
var Ls_inverse = InverseMatrix(Ls);
var r1 = Ls_inverse * yS;
var r2 = ;

function metrics(x1, x2) {
	var s = 0;
	for (var i = 0; i < x1.length; i++) {
		s += Math.pow((x1[i] - x2[i]), 2);
	}
	s = Math.sqrt(s);
	return s;
}


function incasSVM(X, c) {
	var w, w0;
	var ans = [];

	var X1 = [], X2 = [];
	var I_s = [];
	var i1, i2;

	//
	for (var n = 0; n < X.length; n++) {
		if (X[n][1] == 1)
			X1[X1.length] = X[n];
		else
			X2[X2.length] = X[n];
	}

	//
	i1 = X1[0][0];
	var minX1X2 =  metrics(i1, X2[0][0]);
	for (var m = 1; m < X2.length; m++) {
		var buf = metrics(i1, X2[m][0]);
		if (buf < minX1X2) {
			minX1X2 = buf;
			i2 = X2[m][0];
		}
	}
	minX1X2 =  metrics(i2, X1[0][0]);
	for (var m = 1; m < X1.length; m++) {
		var buf = metrics(i2, X1[m][0]);
		if (buf < minX1X2) {
			minX1X2 = buf;
			i1 = X1[m][0];
		}
	}

	//
	

	ans[0] = w;
	ans[1] = w0;
	return ans;
}