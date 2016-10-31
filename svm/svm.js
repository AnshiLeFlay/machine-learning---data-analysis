/*function scalar(u, v) {

}*/
function kernels(u, v, type) {
	if (type == 0) { //<u, v>^2, X=R^2
		var ans = Math.pow(u[0], 2) * Math.pow(v[0], 2) + Math.pow(u[1], 2) * Math.pow(v[1], 2) + 2 * u[0] * v[0] * u[1] * v[1];
	}
	return ans;
}

var Q = [
	[2, -1],
	[-1, 2]
];

function cholesky(Qss) {
	var n = Qss.length;
	var ret = [];
	/*
	double[][] l = new double[m][m];
	for (int i = 0; i< m;i++) {
		for (int k = 0; k < (i+1); k++) {
			double sum = 0;
			for (int j = 0; j < k; j++) {
				sum += l[i][j] * l[k][j];
			}

			l[i][k] = (i == k) ? Math.sqrt(a[i][i] - sum) :
			(1.0 / l[k][k] * (a[i][k] - sum));
		}
	}
	*/
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

console.log(cholesky(Q));