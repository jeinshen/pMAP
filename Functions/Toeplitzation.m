
function Matres=Toeplitzation(Mat)
%this function returns a Toeplitz matrix, closest to Mat for the Frobenius norm.
%this is done by simply averaging along the diagonals.
	[height,width]=size(Mat);  %height>=width required
	Matres=Mat;
	for indexcol=2:width
		valdiag=0;
		valdiag2=0;
		for indexrow=1:width-indexcol+1
			valdiag=valdiag+Mat(indexrow,indexcol-1+indexrow);
			valdiag2=valdiag2+Mat(height-indexrow+1,width-indexcol+2-indexrow);
		end
		valdiag=valdiag/(width-indexcol+1);
		valdiag2=valdiag2/(width-indexcol+1);
		for indexrow=1:width-indexcol+1
			Matres(indexrow,indexcol-1+indexrow)=valdiag;
			Matres(height-indexrow+1,width-indexcol+2-indexrow)=valdiag2;
		end
	end
	for indexcol=1:height-width+1
		valdiag=0;
		for indexrow=1:width
			valdiag=valdiag+Mat(indexcol+indexrow-1,indexrow);
		end
		valdiag=valdiag/width;
		for indexrow=1:width
			Matres(indexcol+indexrow-1,indexrow)=valdiag;
		end
	end
end		%of the function	