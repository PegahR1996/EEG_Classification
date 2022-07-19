%classify for two type condision:delay(decrease speed) and advance(increase speed)
subject=['a','b','c','d','e','i'];
num_sel_f=6;
for iter=8:9;
    num_sample =(iter+1)*10;
    for sbj=1:size(subject,2)
        clear o*
        filename=['C:\Users\10\Documents\MATLAB\mymatlabcode\collectdata\cnt_de_',subject(sbj),'.mat'];
        load(filename);
        filename=['C:\Users\10\Documents\MATLAB\mymatlabcode\collectdata\cnt_ad_',subject(sbj),'.mat'];
        load(filename);
        filename=['C:\Users\10\Documents\MATLAB\mymatlabcode\collectdata\cnt_de_to_prefer_',subject(sbj),'.mat'];
        load(filename);
        filename=['C:\Users\10\Documents\MATLAB\mymatlabcode\collectdata\cnt_ad_to_prefer_',subject(sbj),'.mat'];
        load(filename);
        %ocnt_delay=eeg[onum_channel,(onum_delay*3073)]
        %ocnt_advance=eeg[onum_channel,(onum_advance*3073)]
        cnt_delay=double(cnt_delay);
        cnt_advance=double(cnt_advance);
        cnt_delay_to_prefer=double(cnt_delay_to_prefer);
        cnt_advance_to_prefer=double(cnt_advance_to_prefer);
        %the number of sample in each trial includin 6 secod
        num_s=3073;
        %onum_dealay,onum_advance:the number of trial in delay and advance group
        EEG_constanst=[cnt_delay,cnt_advance];
        EEG_to_prefer=[cnt_delay_to_prefer,cnt_advance_to_prefer];
        number_c=length(EEG_constanst)/num_s;
        number_to_p=length(EEG_to_prefer)/num_s;
        %the number of channel
        num_channel=size(EEG_constanst,1);
        %the number of trial in two group
        %sample rate
        num=number_c + number_to_p;
        fs=512;
        oaccuracy=zeros(1,num);
        
        %make 3d matrix of odata_EEG
        for oj = 0:(number_c -1)
            ojk = num_s*(oj+1);
            data_constant(:,1:num_s,(oj+1))=EEG_constanst(:,oj*num_s+1:ojk);
        end
        
        for oj = 0:(number_to_p -1)
            ojk = num_s*(oj+1);
            data_to_prefer(:,1:num_s,(oj+1))=EEG_to_prefer(:,oj*num_s+1:ojk);
        end
        
        
         d_constant= data_constant(:,end+1-num_sample:end,:);
         d_to_prefer= data_to_prefer(:,102:(num_sample+101),:);
         odata_EEG_3D=cat(3,d_constant,d_to_prefer);
        
        %% leave one out cross validation
        
        for i=1:num
            oind_test = i;
            oind_train = 1:num;
            oind_train(oind_test)=[];
            
            odata_train= zeros(num_channel,num_s,size(oind_train,2));
            odata_test= zeros(num_channel,num_s,1);
            %divide odata_EEG_3D to test and train
            odata_train= odata_EEG_3D(:,:,oind_train);
            odata_test= odata_EEG_3D(:,:,oind_test);
            %making lable
            olable=[ones(1,number_c), 2*ones(1,number_to_p)];
            olable_test=olable(oind_test);
            olable_train =olable(oind_train);
            
            oindex_constant=find(olable_train==1);
            oindex_to_prefer=find(olable_train==2);
            
            num_constant=size(oindex_constant,2);
            num_to_prefer=size(oindex_to_prefer,2);
            
            otrain_constant=zeros(num_channel,num_s,num_constant);
            otrain_to_prefer=zeros(num_channel,num_s,num_to_prefer);
            
            otrain_constant=odata_train(:,:,oindex_constant);
            otrain_to_prefer=odata_train(:,:,oindex_to_prefer);
            %we are not going to use all of 6 second(3073 sample) of data.we
            %are going to use from 52th sample(becouse of 100ms delay in
            %brain commad percieve)to 10 sample after that(52) and add 10 other sample in later
            %iterations of iter
            %100ms/1s * 512 =~52
            
            ofeature_constant=zeros(18,num_constant);
            ofeature_to_prefer=zeros(18,num_to_prefer);
            ofeature_test=zeros(18,1);
            
            ofeature_constant_var=zeros(18,num_constant);
            ofeature_to_prefer_var=zeros(18,num_to_prefer);
            ofeature_test_var=zeros(18,1);
            
            odata_train_constant=zeros(num_channel,num_sample*num_constant);
            odata_train_to_prefer=zeros(num_channel,num_sample*num_to_prefer);
            odata_test_un=zeros(num_channel,num_sample);
            
            %%reshape 3d_eeg data to 2d
            odata_train_constant=reshape(otrain_constant,num_channel,num_sample*(num_constant));
            odata_train_to_prefer=reshape(otrain_to_prefer,num_channel,num_sample*(num_to_prefer));
            odata_test_un=reshape(odata_test,num_channel,num_sample);
            
            
            %% FBCSP
            bands=[4:4:36;8:4:40];
            
            for bn=1:size(bands,2);
                wn=bands(:,bn);
                
                [qq,a]=butter(3,wn/(fs/2),'bandpass');
                
                train_constant_band=zeros(num_channel,num_sample*num_constant);
                train_to_prefer_band=zeros(num_channel,num_sample*num_to_prefer);
                test_band=zeros(num_channel,num_sample);
                
                train_constant_band=filtfilt(qq,a,odata_train_constant);
                train_to_prefer_band=filtfilt(qq,a,odata_train_to_prefer);
                test_band=filtfilt(qq,a,odata_test_un);
                
                %% select channel
                
                feature_de= myvar_all(num_sample,train_constant_band);
                feature_ad= myvar_all(num_sample,train_to_prefer_band);
                
                [H,P]= ttest2(feature_de',feature_ad','Alpha',0.05);
                
                [P,indx]= sort(P,'ascend');
                num_sel_ch=25;
                
                train_constant_band_sel=train_constant_band(indx(1:num_sel_ch),:);
                train_to_prefer_band_sel=train_to_prefer_band(indx(1:num_sel_ch),:);
                test_band_sel= test_band(indx(1:num_sel_ch),:);
                
                
                %calculate the projection matrix for reducing channel
                %number from 30 to 2
                [W] =mycsp(num_sel_ch,num_sample,num_constant,num_to_prefer,train_constant_band_sel,train_to_prefer_band_sel);
                
                w_train_constant =zeros(2,num_sample*num_constant);
                w_train_to_prefer =zeros(2,num_sample*num_to_prefer);
                w_test=zeros(2,num_sample);
                
                w_train_constant =myapplaycsp(num_sample,W,train_constant_band_sel);
                w_train_to_prefer =myapplaycsp(num_sample,W,train_to_prefer_band_sel);
                w_test= myapplaycsp(num_sample,W,test_band_sel);
                %%kamel nashode
                
                f_train_constant_var = zeros(2,num_constant);
                f_train_to_prefer_var =zeros(2,num_to_prefer);
                f_test_var =zeros(2,1);
                
                
                f_train_constant_var = myvar(num_sample,w_train_constant);
                f_train_to_prefer_var =myvar(num_sample,w_train_to_prefer);
                f_test_var =myvar(num_sample,w_test);
                
                q=2*bn;
                qq=(bn-1)*2+1;
                
                ofeature_constant_var(qq:q,:)=f_train_constant_var;
                ofeature_to_prefer_var(qq:q,:)=f_train_to_prefer_var;
                ofeature_test_var(qq:q,:)=f_test_var;
            end
            
            ofeature_constant=ofeature_constant_var;
            ofeature_to_prefer=ofeature_to_prefer_var;
            ofeature_test=ofeature_test_var;
            %% select feature
            
            for fn= 1:18
                om1= mean(ofeature_constant(fn,:));
                om2= mean(ofeature_to_prefer(fn,:));
                os1= var(ofeature_constant(fn,:));
                os2= var(ofeature_to_prefer(fn,:));
                ofdr(fn) = ((om1-om2)^2) / (os1+os2);
            end
            
            [ofdr,oindd]= sort(ofdr,'descend');
            osel_ind = oindd(1:num_sel_f);
            
            ofeature_constant=ofeature_constant(osel_ind,:);
            ofeature_to_prefer=ofeature_to_prefer(osel_ind,:);
            ofeature_test=ofeature_test(osel_ind,:);
            
            ofeature_train= [ofeature_constant,ofeature_to_prefer];
            %         output1=knnclassify(ofeature_test',ofeature_train',olable_train,5);
            omodel= svmtrain(ofeature_train',olable_train,'kernel_function','rbf');
            output1= svmclassify(omodel,ofeature_test');
            %         output1=classify(ofeature_test',ofeature_train',olable_train);
            if olable_test==1
                TP(i)=(olable_test==output1')*100;
                TN(i)=NaN;
            elseif olable_test==2
                TN(i)=(olable_test==output1')*100;
                TP(i)=NaN;
            end
        end
        indP=find(isnan(TP)==1);
        indN=find(isnan(TN)==1);
        TP(indP)=[];
        TN(indN)=[];
        TP_m=mean(TP);
        TN_m=mean(TN);
        accuracy(sbj)=(TP_m*(number_c)+TN_m*(number_to_p))/(number_c+number_to_p);
        sensivity(sbj)=(TP_m*number_c)/(number_c+number_to_p)*2;
        specifity(sbj)=(TN_m*number_to_p)/(number_c+number_to_p)*2;
        clear TP TN
        
        ohh=(num_sample+52)/512;
        disp(['accuracy ',subject(sbj),': ',num2str(accuracy(sbj)),': ',num2str(ohh),':',num2str(num_sel_f),])
        disp(['sensivity ',subject(sbj),': ',num2str(sensivity(sbj)),': ',num2str(ohh),':',num2str(num_sel_f),])
        disp(['specifity ',subject(sbj),': ',num2str(specifity(sbj)),': ',num2str(ohh),':',num2str(num_sel_f),])
        
        
        plot(accuracy(1:sbj),'linewidth',2)
        hold on
        plot(accuracy(1:sbj),'or','linewidth',2)
        grid on
        grid minor
        drawnow
        
    end
    
    disp(['total Accuracy of all subject: ',num2str(mean(accuracy)),':',num2str(ohh),':',num2str(num_sel_f),])
    disp(['total sensivity1 of all subject: ',num2str(mean(sensivity)),':',num2str(ohh),':',num2str(num_sel_f),])
    disp(['total Accuracy of all subject: ',num2str(mean(specifity)),':',num2str(ohh),':',num2str(num_sel_f),])
end

