clear;clc;
tic
path(path,'/home/boonyc/CODES/SegyMAT');

% Required parameters
nt = 3000;  % depend on data traces
ns = 200;
dt = 0.001;
convert_15_to_10 = 0;
convert_15_to_5 = 0;
convert_15_to_2_5 = 1;
wavelet_dir = '/d/tdraid/chaiwoot/Wavelet/DT0.001/';
data_dir = '/d/tdraid/chaiwoot/Modeling/TestModel1/';



[w15,h1,h2] = ReadSu([wavelet_dir 'src15.su']);
[w10,h1,h2] = ReadSu([wavelet_dir 'src10.su']);
[w5,h1,h2] = ReadSu([wavelet_dir 'src5.su']);
[w2_5,h1,h2] = ReadSu([wavelet_dir 'src2.5.su']);
f0 = 15;w15 = w15';
f1 = 10;w10 = w10';
f2 = 5;w5 = w5';
f3 = 2.5;w2_5 = w2_5';

w15f = fft(w15,nt);
w10f = fft(w10,nt);
w5f = fft(w5,nt);
w2_5f = fft(w2_5,nt);

% 10-to-5 Hz Deconvolution
f15_to_10 = w10f.*conj(w15f)./(w15f.*conj(w15f)+1e-10);
f15_to_5 = w5f.*conj(w15f)./(w15f.*conj(w15f)+1e-10);
f15_to_2_5 = w2_5f.*conj(w15f)./(w15f.*conj(w15f)+1e-10);
f5_to_2_5 = w2_5f.*conj(w5f)./(w5f.*conj(w5f)+1e-5);

if convert_15_to_10 == 1
    halfwidth = 200;
    h = blackman(2*halfwidth);
    b1 = h(1:halfwidth);
    b2 = h(halfwidth+1:end);

    for i=1:ns,
        disp(['Shot ' num2str(i)]);
        file_name = [data_dir 'F15_ewt/CSG.'];
        [data,h1,h2] = ReadSu([file_name num2str(i)]);
        ng = size(data,2);
        for ig=1:ng,
            trace = data(:,ig)';
            trace(1:halfwidth) = trace(1:halfwidth).*b1';
            trace(end-halfwidth+1:end) = trace(end-halfwidth+1:end).*b2';
%            plot(trace);
            trace = real(ifft(fft(trace).*f15_to_10));
            data(:,ig) = trace';
%            hold on;plot(trace,'r');hold off;pause;
        end;
        file_name = [data_dir 'F10/CSG.'];
        WriteSu([file_name num2str(i)],data,'dt',h1(1).dt/1e6);
    end;
end

if convert_15_to_5 == 1
    halfwidth = 200;
    h = blackman(2*halfwidth);
    b1 = h(1:halfwidth);
    b2 = h(halfwidth+1:end);

    for i=1:ns,
        disp(['Shot ' num2str(i)]);
        file_name = [data_dir 'F15_ewt/CSG.'];
        [data,h1,h2] = ReadSu([file_name num2str(i)]);
        ng = size(data,2);
        for ig=1:ng,
            trace = data(:,ig)';
            trace(1:halfwidth) = trace(1:halfwidth).*b1';
            trace(end-halfwidth+1:end) = trace(end-halfwidth+1:end).*b2';
%            plot(trace);
            trace = real(ifft(fft(trace,nt).*f15_to_5));
            trace(1:halfwidth) = trace(1:halfwidth).*b1';
            trace(end-halfwidth+1:end) = trace(end-halfwidth+1:end).*b2';
            data(:,ig) = trace';
%            hold on;plot(data(:,ig),'r');hold off;pause;
        end;
        file_name = [data_dir 'F5/CSG.'];
        WriteSu([file_name num2str(i)],data,'dt',h1(1).dt/1e6);
    end;
end

if convert_15_to_2_5 == 1
    halfwidth = 400;
    h = blackman(2*halfwidth);
    b1 = h(1:halfwidth);
    b2 = h(halfwidth+1:end);

    for i=1:ns,
        disp(['Shot ' num2str(i)]);
        file_name = [data_dir 'F15_ewt/CSG.'];
        [data,h1,h2] = ReadSu([file_name num2str(i)]);
        ng = size(data,2);
        for ig=1:ng,
            trace = data(:,ig)';
            trace(1:halfwidth) = trace(1:halfwidth).*b1';
            trace(end-halfwidth+1:end) = trace(end-halfwidth+1:end).*b2';
%            plot(trace);
            trace = real(ifft(fft(trace,nt).*f15_to_2_5));
            trace(1:halfwidth) = trace(1:halfwidth).*b1';
            trace(end-halfwidth+1:end) = trace(end-halfwidth+1:end).*b2';
            data(:,ig) = trace';
%            hold on;plot(data(:,ig),'r');hold off;pause;
        end;
        file_name = [data_dir 'F2.5/CSG.'];
        WriteSu([file_name num2str(i)],data,'dt',h1(1).dt/1e6);
    end;
end

disp(['Runtime = ' num2str(toc/60) ' min']);
