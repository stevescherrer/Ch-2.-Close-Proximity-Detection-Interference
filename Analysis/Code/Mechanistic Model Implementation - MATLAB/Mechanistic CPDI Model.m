%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Name: Brendan Rideout (University of Hawaii - Ocean and Resources Engineering)
%
%   Course: CPDI Modeling
%
%   Date: 25/06/2017
%
%   Description: Linear trajectory acoustic propagation forward model
%
%   Inputs:		m               (vector of model parameter values)
%				arrivePath      (binary vector indicating which arrival paths arrival times must be calculated for)
%
%   Outputs:	tArrivalMod     (model predicted arrival times)
%               lTravel         (propagation length for each arrival)
%              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tArrivalMod, lTravel] = MOI_XYZ(m, arrivePath)
 
    xSource = m(1);
    ySource = m(2);
	zSource = m(3);
    cBartSource = m(4);
    cBartOffset = m(5);
	xPhone = m(6);
    yPhone = m(7);
    zPhone = m(8);
	dWater = m(9);
	cWater = m(10);
    if(length(m)>=11)
        cConstant = m(11);
        cConstantFLAG = true;
    end
	
	D = arrivePath(1);
	B = arrivePath(2);
	T = arrivePath(3);
	TB = arrivePath(4);
	
    BT = arrivePath(5);
    BTB = arrivePath(6);
    TBT = arrivePath(7);
    TBTB = arrivePath(8);
    
    BTBT = arrivePath(9);
    BTBTB = arrivePath(10);
    TBTBT = arrivePath(11);
    TBTBTB = arrivePath(12);
	
	BTBTBT = arrivePath(13);
	BTBTBTB = arrivePath(14);
	TBTBTBT = arrivePath(15);
	TBTBTBTB = arrivePath(16);
    
    BTBTBTBT = arrivePath(17);
    BTBTBTBTB = arrivePath(18);
    TBTBTBTBT = arrivePath(19);
    TBTBTBTBTB = arrivePath(20);
	
	delX = (xSource - xPhone)^2;
    delY = (ySource - yPhone)^2;
    delR = delX + delY;
	
	ndat = sum(sum(arrivePath)); %single hydrophone case
	lTravel = zeros(ndat,1);
	currIndex = 0;
	
	if(D) %direct
		currIndex = currIndex + 1;
		nBottBounc = 0;
		lTravel(currIndex,1) = sqrt(delR + (2*nBottBounc*dWater + zPhone - zSource)^2);
	end
	
		if(B) %bottom
			currIndex = currIndex + 1;
			nBottBounc = 1;
			lTravel(currIndex,1) = sqrt(delR + (2*nBottBounc*dWater - zPhone - zSource)^2);
		end
		
		if(T) %top
			currIndex = currIndex + 1;
			nBottBounc = 0;
			lTravel(currIndex,1) = sqrt(delR + (2*nBottBounc*dWater + zPhone + zSource)^2);
		end
		
		if(TB) %topBott
			currIndex = currIndex + 1;
			nBottBounc = 1;
			lTravel(currIndex,1) = sqrt(delR + (2*nBottBounc*dWater - zPhone + zSource)^2);
		end
	
	if(BT) %bottTop
		currIndex = currIndex + 1;
		nBottBounc = 1;
		lTravel(currIndex,1) = sqrt(delR + (2*nBottBounc*dWater + zPhone - zSource)^2);
    end
    
		if(BTB) %bottTopBott
			currIndex = currIndex + 1;
			nBottBounc = 2;
			lTravel(currIndex,1) = sqrt(delR + (2*nBottBounc*dWater - zSource - zPhone)^2);
		end        
		
		if(TBT) %topBottTop
			currIndex = currIndex + 1;
			nBottBounc = 1;
			lTravel(currIndex,1) = sqrt(delR + (2*nBottBounc*dWater + zPhone + zSource)^2);
		end
		
		if(TBTB) %topBottTopBott
			currIndex = currIndex + 1;
			nBottBounc = 2;
			lTravel(currIndex,1) = sqrt(delR + (2*nBottBounc*dWater - zPhone + zSource)^2);
		end
    	
	if(BTBT) %bottTopBottTop
		currIndex = currIndex + 1;
		nBottBounc = 2;
		lTravel(currIndex,1) = sqrt(delR + (2*nBottBounc*dWater + zPhone - zSource)^2);
    end
    
		if(BTBTB) %bottTopBottTopBott
			currIndex = currIndex + 1;
			nBottBounc = 3;
			lTravel(currIndex,1) = sqrt(delR + (2*nBottBounc*dWater - zSource - zPhone)^2);
		end
		
		if(TBTBT) %topBottTopBottTop
			currIndex = currIndex + 1;
			nBottBounc = 2;
			lTravel(currIndex,1) = sqrt(delR + (2*nBottBounc*dWater + zPhone + zSource)^2);
		end
		
		if(TBTBTB) %topBottTopBottTopBott
			currIndex = currIndex + 1;
			nBottBounc = 3;
			lTravel(currIndex,1) = sqrt(delR + (2*nBottBounc*dWater - zPhone + zSource)^2);
		end
	
	if(BTBTBT) %bottTopBottTopBottTop
		currIndex = currIndex + 1;
		nBottBounc = 3;
		lTravel(currIndex,1) = sqrt(delR + (2*nBottBounc*dWater + zPhone - zSource)^2);
    end
	
        if(BTBTBTB)
            currIndex = currIndex + 1;
            nBottBounc = 4;
            lTravel(currIndex,1) = sqrt(delR + (2*nBottBounc*dWater - zSource - zPhone)^2);
        end
    
		if(TBTBTBT) %topBottTopBottTopBottTop
			currIndex = currIndex + 1;
			nBottBounc = 3;
			lTravel(currIndex,1) = sqrt(delR + (2*nBottBounc*dWater + zPhone + zSource)^2);
        end
        
        if(TBTBTBTB)
            currIndex = currIndex + 1;
			nBottBounc = 4;
			lTravel(currIndex,1) = sqrt(delR + (2*nBottBounc*dWater - zPhone + zSource)^2);
        end
		
	if(BTBTBTBT)
		currIndex = currIndex + 1;
		nBottBounc = 4;
		lTravel(currIndex,1) = sqrt(delR + (2*nBottBounc*dWater + zPhone - zSource)^2);
    end
        
        if(BTBTBTBTB)
            currIndex = currIndex + 1;
            nBottBounc = 5;
            lTravel(currIndex,1) = sqrt(delR + (2*nBottBounc*dWater - zSource - zPhone)^2);
        end
	
		if(TBTBTBTBT)
			currIndex = currIndex + 1;
			nBottBounc = 4;
			lTravel(currIndex,1) = sqrt(delR + (2*nBottBounc*dWater + zPhone + zSource)^2);
        end
        
        if(TBTBTBTBTB)
            currIndex = currIndex + 1;
			nBottBounc = 5;
			lTravel(currIndex,1) = sqrt(delR + (2*nBottBounc*dWater - zPhone + zSource)^2);
        end                                  
	
	
	tTravelMod = lTravel / cWater;
    if(cConstantFLAG)
        tArrivalMod = tTravelMod + (1/cConstant)*(cBartSource) + (1/cConstant)*(cBartOffset);
    else
        tArrivalMod = tTravelMod + tSource + tOffset;
    end