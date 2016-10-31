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

	return ans;
}

//solving the quadratic subproblem
//input data 
//X_l, C
//Q_ss, y_s, e, c_c
var Xl = [];
var C = 1;



//2 closest points
//for example
var I_s = [];
var I_o = [];
var I_c = [];