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
	//need optimization this
	for (var i = 0; i < ret.length; i++) {
		for (var j = 0; j < ret.length; j++) {
			if (typeof(ret[i][j]) == 'undefined') ret[i][j] = 0;
		}
	}

	return ret;
}

function transp(mtrx) {
	var ans = [];
	var m = mtrx.length;
	var n = mtrx[0].length;
	var vector = false;
	if (typeof(n) == 'undefined') {
		vector = true;
		n = 1;
	}
	for (var i = 0; i < n; i++) {
		ans[i] = [];
		for (var j = 0; j < m; j++) {
			if (vector == true) {
				ans[j] = [];
				ans[j][i] = mtrx[j];
			}
			else ans[i][j] = mtrx[j][i];
		}
	}
	return ans;
}

function mult(mtrxA, mtrxB) {
	var ans = [];

	var l = mtrxA.length;
	var m = mtrxA[0].length;

	var mb = mtrxB.length;
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
function plusMtrx(mtrxA, mtrxB, sign) {
	var ans = [];
	if (mtrxA[1].length > 1) {
		for (var i = 0; i < mtrxA.length; i++) {
			ans[i] = [];
			for (var j = 0; j < mtrxA[0].length; j++) {
				if (sign == '+') {
					if (typeof(mtrxB) == 'number')
						ans[i][j] = mtrxA[i][j] + mtrxB;
					else
						ans[i][j] = mtrxA[i][j] + mtrxB[i][j];
				} else {
					if (typeof(mtrxB) == 'number')
						ans[i][j] = mtrxA[i][j] - mtrxB;
					else
						ans[i][j] = mtrxA[i][j] - mtrxB[i][j];
				}
			}
		}
	} else {
		for (var i = 0; i < mtrxA.length; i++) {
			ans[i] = [];
			if (sign == '+') {
				if (typeof(mtrxB) == 'number')
					ans[i][0] = mtrxA[i] + mtrxB;
				else
					ans[i][0] = mtrxA[i] + mtrxB[i];
			} else {
				if (typeof(mtrxB) == 'number')
					ans[i][0] = mtrxA[i] - mtrxB;
				else
					ans[i][0] = mtrxA[i] - mtrxB[i];
			}
		}
	}
	return ans;
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
    //if (typeof(A[0]) == 'undefined') return A[0];

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
	if (A.length == 1 && A[0].length == 1) 
	{
		var ret = [];
		ret[0] = []; 
		ret[0][0] = 1 / A[0];
		return ret;
	} 
    var det = Determinant(A);
    if (det == 0) return false;
    var N = A.length, A = AdjugateMatrix(A);
    for (var i = 0; i < N; i++)
     { for (var j = 0; j < N; j++) A[i][j] /= det; }
    return A;
}
//////////////////////////////////////////////////////////////////////////////////////////////

function mtrxLog(mtrx, msg) {
	console.log(msg);
	for (var i = 0; i < mtrx.length; i++)
		console.log(mtrx[i]);

}

function metrics(x1, x2) {
	var s = 0;
	for (var i = 0; i < x1.length; i++) {
		s += Math.pow((x1[i] - x2[i]), 2);
	}
	s = Math.sqrt(s);
	return s;
}


function incasSVM(X, c) {
	var w = [0, 0], w0 = 0;
	var ans = [];

	var X1 = [], X2 = [];
	var I_s = [], I_o = [], I_c = [];
	var i1, i2;

	var flag = 0, flag2 = 0;

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
	var del = 0;

	for (var m = 1; m < X2.length; m++) {
		var buf = metrics(i1, X2[m][0]);
		if (buf < minX1X2) {
			minX1X2 = buf;
			i2 = X2[m][0];
			del = m;
		}
	}
	X2.splice(del, 1);
	minX1X2 =  metrics(i2, X1[0][0]);
	del = 0;
	for (var m = 1; m < X1.length; m++) {
		var buf = metrics(i2, X1[m][0]);
		if (buf < minX1X2) {
			minX1X2 = buf;
			i1 = X1[m][0];
			del = m;
		}
	}
	X1.splice(del, 1);
	//
	
	I_s[I_s.length] = [i1, 1];
	I_s[I_s.length] = [i2, -1];

	I_o = I_o.concat(X1);
	I_o = I_o.concat(X2);
	do {
		flag2 = 0;
		do {
			flag = 0;
			var Qss = [];
			var Qcs = [];
			//init Qss, Qcs
			var modI_s = I_s.length;
			var modI_c = I_c.length;
			for (var i = 0; i < modI_s; i++) {
				Qss[i] = [];
				for (var j = 0; j < modI_s; j++) {
					Qss[i][j] = I_s[i][1] * I_s[j][1] * kernels(I_s[i][0], I_s[j][0], 0);
				}
			}
			mtrxLog(Qss, 'Qss test');
			var L_s = cholesky(Qss);
			mtrxLog(L_s, 'L_s test');

			var yC = [];
			var eS = [];

			//init yS, eS
			var yS = [];
			var eS = [];
			for (var i = 0; i < modI_s; i++) {
				yS[i] = I_s[i][1];
				eS[i] = -1 * c;
			}

			yS = transp(yS);
			eS = transp(eS);

			//
			var eCC = [];
			for (var i = 0; i < modI_c; i++) {
				eCC[i] = c;
			}

			//solving quadratic porblem
			var Ls_inverse = InverseMatrix(L_s);
			var r1 = mult(Ls_inverse, yS);
			var r2;
			var beta;
			var alphaS;
			if (I_c.length != 0) {
				r2 = mult(transp(eCC), Qcs);
				r2 = transp(plusMtrx(eS, r2));
				r2 = mult(Ls_inverse, r2);

				beta = plusMtrx(mult(transp(r1), r2), mult(tranps(eCC), yC), '-');
				beta = mult(beta, InverseMatrix(mult(transp(r1), r1)));
			} else {
				r2 = mult(Ls_inverse, eS);

				beta = mult(mult(transp(r1), r2), InverseMatrix(mult(transp(r1), r1)));
			}

			alphaS = mult(transp(Ls_inverse), plusMtrx(mult(r1, beta), r2, '-'));

			//ifs 
			//if flag == 1 go to solving qp
			//if 
			if (modI_s > 2) {
				for (var i = 0; i < modI_s; i++) {
					if (alphaS[i][0] <= 0) {
						I_o.push(I_s[i]);
						I_s.splice(i, 1);
						modI_s--;
						flag = 1;
						break;
					}
				}
			}

			if (modI_s > 2 && flag != 1) {
				for (var i = 0; i < modI_s; i++) {
					if (alphaS[i][0] >= c) {
						I_c.push(I_s[i]);
						I_s.splice(i, 1);
						modI_s--;
						flag = 2;
						break;
					}
				}
			}
			console.log(flag);
		} while(flag != 0);
		/*
		do {
			do {

			}
			while();
		}
		while();
		*/

		mtrxLog(I_s, 'I_s : '); //alphaS
		mtrxLog(I_c, 'I_c : '); //C*eC
		mtrxLog(I_o, 'I_o : ');	//0
		mtrxLog(alphaS, 'alphaS : ');

		w[0] = 0;
		w[1] = 0;
		w0 = 0;


		for (var i = 0; i < I_s.length; i++) {
			w[0] += alphaS[i][0] * I_s[i][1] * I_s[i][0][0];
			w[1] += alphaS[i][0] * I_s[i][1] * I_s[i][0][0];
		}
		for (var i = 0; i < I_c.length; i++) {
			w[0] += c * I_c[i][1] * I_c[i][0][0];
			w[1] += c * I_c[i][1] * I_c[i][0][0];
		}


		var sca = 0;
		w0 = 0;
		for (var i = 0; i < X.length; i++) {
			sca = 0;
			for (var j = 0; j < X[i][0].length; j++) {
				sca += w[j] * X[i][0][j];
				console.log(w[j]);
			}
			w0 += sca - X[i][1];
		}
		w0 = w0 / X.length;
		console.log('w0');
		console.log(w0);

		var marginI_o = [];
		var marginI_c = [];
		
		for (var i = 0; i < I_o.length; i++) {
			sca = 0;
			for (var j = 0; j < I_o[i][0].length; j++) {
				sca += w[j] * I_o[i][0][j];
			} 
			marginI_o.push(sca - w0);
		}

		for (var i = 0; i < I_c.length; i++) {
			sca = 0;
			for (var j = 0; j < I_c[i][0].length; j++) {
				sca += w[j] * I_c[i][0][j];
			} 
			marginI_c.push(sca - w0);
		}
		console.log('I_o 1');
		console.log(I_o);
		console.log('I_c 1');
		console.log(I_c);

		for (var i = 0; i < marginI_o.length; i++) {
			if (marginI_o[i] <= 1) {
				I_s.push(I_o[i]);
				I_o.splice(i, 1);
				flag2++;
				break;
			}
		}

		for (var i = 0; i < marginI_c.length; i++) {
			if (marginI_c[i] >= 1) {
				I_s.push(I_c[i]);
				I_c.splice(i, 1);
				flag2++;
				break;
			}
		}

		console.log('I_o 2');
		console.log(I_o);
		console.log('I_c 2');
		console.log(I_c);

		console.log(flag2);
	} while(flag2 != 0);

	ans[0] = w;
	ans[1] = w0;

	return ans;
}