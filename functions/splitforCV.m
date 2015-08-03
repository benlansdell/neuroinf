function [data, testdata] = splitforCV(data, percentage)
    nB = size(data.X, 1);
    trainidx = ceil(nB*percentage/100);
    testdata = data;
    data.X = data.X(1:trainidx,:);
    data.y = data.y(1,1:trainidx);
	data.cursor = data.cursor(1:trainidx,:); 
	data.grip = data.grip(1:trainidx,:);
    testdata.X = testdata.X(trainidx+1:end,:);
    testdata.y = testdata.y(trainidx+1:end,:);
	testdata.cursor = testdata.cursor(trainidx+1:end,:); 
	testdata.grip = testdata.grip(trainidx+1:end,:);
end