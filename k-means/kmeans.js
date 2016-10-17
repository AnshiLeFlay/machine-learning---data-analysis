//0 - Euclidean distance
function metrics(x1, x2, type) {
	var s = 0;
	if (x1.length == x2.length) {
		if (type == 0) {
			for (var i = 0; i < x1.length; i++) {
				s += Math.pow((x1[i] - x2[i]), 2);
			}
			s = Math.sqrt(s);
			return s;
		} else {
			s = 'err2';
			return s;
		}
	} else {
		s = 'err1';
		return s;
	} 
}

function minArrIndexes(arr) {
	var indexes = [];
	var arr_min = Math.min.apply(null, arr);
	var f1 = arr.indexOf(arr_min);
	indexes[0] = f1;
	for (var i = 0; i < arr.length; i++) {
		if (arr.indexOf(arr_min, parseInt(f1)+1) !== -1) {
			f1 = arr.indexOf(arr_min, parseInt(f1)+1);
			indexes[indexes.length] = f1;
		} else {
			break;	
		}	
	}
	return indexes;
}

function randomInt(min, max) {
	var rand = min - 0.5 + Math.random() * (max - min + 1)
	rand = Math.round(rand);
	return rand;
}

function kMeans(vectorM, k) {
	//количество кластеров

	//данные для обработки

	//выбор начальных центров
	var centers = [];

	for (var i = 0; i < k; i++) {
		centers[i] = vectorM[i];
	}
	//console.log("centers : " + centers);
	//console.log(centers[1]);

	//обработка
	for (var steps = 0; steps < 1000; steps++)
	{
		var metrics_buffer = [];
		var mins = [];
		var clusters = new Array(k);

		for (var h = 0; h < k; h++)
			clusters[h] = [];

		for (var m = 0; m < vectorM.length; m++) {
			for (var j = 0; j < k; j++) {
				metrics_buffer[j] = metrics(vectorM[m], centers[j], 0);
			}
			//console.log('distance ' + vectorM[m] + ' to clusters : ' + metrics_buffer);
			mins = minArrIndexes(metrics_buffer);
			//console.log('mins : ' + mins);
			var toCluster = mins[parseInt(randomInt(0, mins.length-1))];
			var l = clusters[parseInt(toCluster)].length;
			//console.log(l);
			clusters[parseInt(toCluster)][l] = vectorM[m]; 
		}

		//console.log('result of 1 step : ');
		//console.log(clusters);
		var x;
		var y;
		for (var i = 0; i < clusters.length; i++) {
			x = 0;
			y = 0;
			for (var j = 0; j < clusters[i].length; j++) {
				x += parseInt(clusters[i][j][0]);
				y += parseInt(clusters[i][j][1]);
			}
			x = x / clusters[i].length;
			y = y / clusters[i].length;
			centers[i] = [x, y];
		}
	}

	//console.log(clusters[][][]);
	return clusters;
}