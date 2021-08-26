function [ ] = EMmod2CSEMResp(varargin)
% convert the data file of EMmod program by Hunziker et al. (2015)
% to CSEMResp_1.0 
%                                         --- HB, Sep. 20, 2018

close all; fclose all; clear; clc;

datInfo.freqs = [0.25];
datInfo.txLoc = [0 0 950 90 0; 0 0 950 0 0];
datInfo.dataType = {'Ex','Ey','Ez','Bx','By','Bz'};

datInfo.dpLen = 0;

% receiver location range
rxRange.x = [0 0];
rxRange.y = [-2000 8000];

nRx = 101;

nFreq = length(datInfo.freqs);
nTx   = size(datInfo.txLoc,1);
nDT   = length(datInfo.dataType);


% initialize the fields of ddata
ddata.rxLoc=[]; 
% ddata.Ex=[];  ddata.Ey=[];  ddata.Ez=[];
% ddata.Bx=[];  ddata.By=[];  ddata.Bz=[];
ddata.Ex=zeros(nRx,nTx,nFreq);
ddata.Ey=zeros(nRx,nTx,nFreq);
ddata.Ez=zeros(nRx,nTx,nFreq);
ddata.Bx=zeros(nRx,nTx,nFreq);
ddata.By=zeros(nRx,nTx,nFreq);
ddata.Bz=zeros(nRx,nTx,nFreq);

% added on 2020-05-08
geometry = cell(nTx,1);
for iTx=1:nTx
    if datInfo.txLoc(iTx, 4) == 90
        geometry{iTx} = 'inline';
    elseif datInfo.txLoc(iTx, 4) == 0
        geometry{iTx} = 'broadside';
    end
end

for iFreq=1:nFreq  
    for iTx=1:nTx
        for iDT=1:nDT
            dt = datInfo.dataType{iDT};
            
            if strcmp(geometry{iTx}, 'inline')  && (strcmp(dt, 'Ex') || strcmp(dt, 'By') || strcmp(dt, 'Bz'))
                continue;
            end
            
            if strcmp(geometry{iTx}, 'broadside')  && (strcmp(dt, 'Ey') || strcmp(dt, 'Ez') || strcmp(dt, 'Bx'))
                continue;
            end
            
            RSLTSpec = {'*.rslt','RSLT Files(*.rslt)'; '*.*', 'All Files(*.*)'};
            [dfile, dpath] = uigetfile(RSLTSpec, ['Open EMmod result file: ',...
                num2str(datInfo.freqs(iFreq)),'Hz,',...
                geometry{iTx},', ',dt]);
            
            if ~ischar(dfile), continue; end

            emmodfile = fullfile(dpath, dfile);
            [dataTmp] = readEMmodFile(emmodfile,rxRange);
            
            if iFreq==1 && iTx==1
                ddata.rxLoc = dataTmp.rxLoc;
            end
            
            if strcmp(dt, 'Ex')
                ddata.Ex(:,iTx,iFreq) = dataTmp.data;
            elseif strcmp(dt, 'Ey')
                ddata.Ey(:,iTx,iFreq) = dataTmp.data;
            elseif strcmp(dt, 'Ez')
                ddata.Ez(:,iTx,iFreq) = dataTmp.data;
            elseif strcmp(dt, 'Bx')
                ddata.Bx(:,iTx,iFreq) = dataTmp.data;
            elseif strcmp(dt, 'By')
                ddata.By(:,iTx,iFreq) = dataTmp.data;
            elseif strcmp(dt, 'Bz')
                ddata.Bz(:,iTx,iFreq) = dataTmp.data;
            end
            
            clear dataTmp;
        end
    end
end

% nRx = size(ddata.rxLoc,1);
% if isempty(ddata.Ex), ddata.Ex=zeros(nRx,nTx,nFreq); end
% if isempty(ddata.Ey), ddata.Ey=zeros(nRx,nTx,nFreq); end
% if isempty(ddata.Ez), ddata.Ez=zeros(nRx,nTx,nFreq); end
% if isempty(ddata.Bx), ddata.Bx=zeros(nRx,nTx,nFreq); end
% if isempty(ddata.By), ddata.By=zeros(nRx,nTx,nFreq); end
% if isempty(ddata.Bz), ddata.Bz=zeros(nRx,nTx,nFreq); end

% from H to B
MU0 = 4*pi*1e-7;
ddata.Bx = MU0 * ddata.Bx;
ddata.By = MU0 * ddata.By;
ddata.Bz = MU0 * ddata.Bz;

datInfo.rxLoc = ddata.rxLoc;

fwdResp = horzcat(ddata.Ex,ddata.Ey,ddata.Ez,ddata.Bx,ddata.By,ddata.Bz);
clear ddata;

[dpath,fname,ext] = fileparts(emmodfile);
[fname,~] = strtok(fname,'_');
outfile = fullfile(dpath,[fname,'_emmod.resp']);
writeCSEMResp(outfile, datInfo, fwdResp);

return;
 end


%==========================================================================
%=========================================================== readEMmodFile 
%==========================================================================
function [data] = readEMmodFile(datafile,rxRange)
% read the data file of EMmod program by Hunziker et al. (2015).
%                                --- HB, June 23, 2016

fid = fopen(datafile, 'r');
fgetl(fid);

j = 0;
while ~feof(fid)
    j = j+1;
    rtmp = fscanf(fid,'%f\n',5);
    mtmp(j,:) = rtmp;
end

x1 = rxRange.x(1);
x2 = rxRange.x(2);
y1 = rxRange.y(1);
y2 = rxRange.y(2);

id = find(mtmp(:,1)>=x1 & mtmp(:,1)<=x2 & mtmp(:,2)>=y1 & mtmp(:,2)<=y2);

data.rxLoc = mtmp(id,1:3);
data.data  = complex(mtmp(id,4), mtmp(id,5));

fclose(fid);
return;
end


%==========================================================================
%=========================================================== writeCSEMResp 
%==========================================================================
function writeCSEMResp(respfile, datInfo, fwdResp)

    fid = fopen(respfile, "wt");

    fprintf(fid,"%-18s %s\n","# Format:","CSEMResp_1.0");
    descrb  = ['Response file generated at ', datestr(now)]; 
    fprintf(fid,"%-18s %s\n","# Description:", descrb);

    % source type
    srcType = "HED";
    fprintf(fid,"%-20s %s\n","Source Type:",srcType);

    % dipole length
    if datInfo.dpLen > 1e-1
        fprintf(fid,"%-20s %G\n","Dipole Length:", datInfo.dpLen);
    end

    % source location
    txLoc = datInfo.txLoc;
    nTx = size(txLoc, 1);
    fprintf(fid,"%-25s %4d\n","Source Location (m):", nTx);
    fprintf(fid,"%-10s %-12s %-12s %-6s %-16s %s\n","#","X","Y","Z","Azimuth","Dip");
    for i = 1:nTx
        for j = 1:5
            fprintf(fid,"%12.2f ",txLoc(i,j));
        end
        fprintf(fid,"\n");
    end

    % receiver location
    rxLoc = datInfo.rxLoc;
    nRx = size(rxLoc, 1);
    fprintf(fid,"%-25s %4d\n","Receiver Location (m):", nRx);
    fprintf(fid,"%-10s %-12s %-12s %s\n","#","X","Y","Z");
    for i = 1:nRx
        fprintf(fid,"%12.2f %12.2f %12.2f\n", rxLoc(i,1), rxLoc(i,2), rxLoc(i,3));
    end

    % frequencies
    freqs = datInfo.freqs;
    nFreq = length(freqs);
    fprintf(fid,"%-20s %4d\n","Frequencies (Hz):", nFreq);
    for i = 1:nFreq
        fprintf(fid,"%15.5e\n", freqs(i));
    end

    % data type
    dataType = datInfo.dataType;
    nDT = length(dataType);
    fprintf(fid,"DataType:  %4d\n", nDT);
    for i = 1:nDT
        fprintf(fid,"%s\n",dataType{i});
    end

    nData  = nFreq * nTx * nRx;

    fprintf(fid,"%-15s %d\n","Data Block:",nData);


    fprintf(fid,"# %-8s %-7s %-17s %-29s %-29s %-29s %-29s %-29s %s\n","FreqNo.","TxNo.","RxNo.",...
            "Ex(Re,Im)","Ey(Re,Im)","Ez(Re,Im)","Bx(Re,Im)","By(Re,Im)","Bz(Re,Im)");

    for i = 1:nFreq
        for j=1:nTx
            for k=1:nRx
                fprintf(fid,"%6d %8d %7d %17.6e %13.6e %15.6e %13.6e %15.6e %13.6e %15.6e %13.6e %15.6e %13.6e %15.6e %13.6e\n",...
                i, j, k, real(fwdResp(k,j,i)), imag(fwdResp(k,j,i)), real(fwdResp(k,nTx+j,i)), imag(fwdResp(k,nTx+j,i)),...
                real(fwdResp(k,2*nTx+j,i)), imag(fwdResp(k,2*nTx+j,i)), real(fwdResp(k,3*nTx+j,i)), imag(fwdResp(k,3*nTx+j,i)),...
                real(fwdResp(k,4*nTx+j,i)), imag(fwdResp(k,4*nTx+j,i)), real(fwdResp(k,5*nTx+j,i)), imag(fwdResp(k,5*nTx+j,i)) );
            end
        end
    end

    fclose(fid);

end