% LOAD FLOW SOLUTION USING THE GAUSS SEIDEL LOAD FLOW METHOD OF NUMERICAL SOLUTION OF
% NON-LINEAR EQUATION(s)

% ASSIGNMENT 2

% LOAD FLOW SOLUTION USING THE FAST-DECOUPLED LOAD FLOW METHOD OF NUMERICAL USING ANGLE IN RADIANT SOLUTION OF

MVA_base=100.0; %Defining the Base-MVA
%Calling the function containing the input data
[linedata,busdata,tol,acceleration_factor,nPV,nLB,lPV] = Bus_line_inputs();


tolerence=tol;
fb=linedata(:,1); 
tb=linedata(:,2); 
r=linedata(:,3); 
x=linedata(:,4); 
b=linedata(:,5); 
a=linedata(:,6); 
bus_type=(busdata(:,2)); 
v_bus_initial=(busdata(:,3)); 
P_gen=((busdata(:,5))*(MVA_base)); 
Q_gen=((busdata(:,6))*(MVA_base)); 
S_gen=((P_gen)+(1i)*(Q_gen)); 
P_load=((busdata(:,7))*(MVA_base));
Q_load=((busdata(:,8))*(MVA_base)); 
S_load=((P_load)+(1i)*(Q_load)); 
Q_min=((busdata(:,9))*(MVA_base)); 
Q_max=((busdata(:,10))*(MVA_base)); 
B_shunt=((busdata(:,11)));
z_elementary=((r)+(1i)*(x));
b_elementary=((1i)*(b)); 
nbus = max(max(fb),max(tb));   
nbranch = length(fb);          

z=zeros(nbus,nbus); 
b=zeros(nbus,nbus); 
Ybus=zeros(nbus,nbus);
tolerence_checker=1; 
iteration=0; 
% Initialising Values
V_new=v_bus_initial;
v_scheduled=v_bus_initial; 
V_new_accelerated=zeros(1,nbus); 
V_new_del=zeros(1,nbus); 
difference=zeros(1,nbus); 
real_diff=zeros(1,nbus); 
imag_diff=zeros(1,nbus); 
delta=zeros(1,nbus); 
E_new=zeros(1,nbus);
F_new=zeros(1,nbus); 
Q_intermediate=zeros(1,nbus); 
Q_final=zeros(1,nbus);
complex_flow_line=zeros(1,nbranch);
line_flows=zeros(nbus,nbus);
active_flow_line=zeros(1,nbranch); 
reactive_flow_line=zeros(1,nbranch); 
bus_power_injection=zeros(1,nbus);
bus_power_mismatch=zeros(1,nbus);
line_loss=zeros(nbus,nbus);
flow_count=1; 
sum_line_loss=((0.0)+(1i)*(0.0)); 
sum=zeros(1,nbus); 
shunt_fb_onr=zeros(1,nbranch); 
shunt_tb_onr=zeros(1,nbranch); 
P_injected_sum=zeros(nbus,1); 
Q_injected_sum=zeros(nbus,1); 
P_injected_bus=zeros(nbus,1); 
Q_injected_bus=zeros(nbus,1); 
partial_P_delta=zeros(nbus-1,nbus-1); 
partial_Q_delta=zeros(nbus-nPV-1,nbus-1); 
partial_P_vol_mag=zeros(nbus-1,nbus-nPV-1); 
partial_Q_vol_mag=zeros(nbus-nPV-1,nbus-nPV-1); 
J=zeros(((2*(nbus-1))-nPV),((2*(nbus-1))-nPV)); 
inj_active_pow_mismatch_vector=zeros(nbus-1,1); 
inj_reactive_pow_mismatch_vector=zeros(nbus-nPV-1,1); 
inj_pow_mismatch_vector=zeros(((2*(nbus-1))-nPV),1); 
correction_vector=zeros(((2*(nbus-1))-nPV),1);
correction_voltage_angle=zeros((nbus-1),1); 
correction_voltage_magnitude=zeros((nbus-1-nPV),1); 
mismatch_active_vector_element_count=1; 
mismatch_reactive_vector_element_count=1; 
mismatch_vector_element_count=0;
real_mismatch=zeros(nbus-1,1); 
imag_mismatch=zeros(nbus-1-nPV,1); 
V_new_for_Q=zeros((nbus),1);
BforP=zeros((nbus-1),(nbus-1));
BforQ=zeros((nbus-1-nPV),(nbus-1-nPV));
I_count_for_BQ=1;
K_count_for_BQ=1;

for u=1:((2*(nbus-1))-nPV)
    correction_vector(u)=1.0; 
end


for u=1:nbranch
    z(fb(u),tb(u))=z_elementary(u); 
    z(tb(u),fb(u))=z(fb(u),tb(u)); 
    b(fb(u),tb(u))=b_elementary(u); 
    b(tb(u),fb(u))=b(fb(u),tb(u));
end

for u=1:nbranch
    if((a(u,1))~=1.0)
        shunt_fb_onr(1,u)=(((a(u,1))^(2))/(1-(a(u,1)))); 
        b(fb(u),tb(u))=(1/((shunt_fb_onr(1,u))*(z(fb(u),tb(u)))))+(b(fb(u),tb(u))); 
        shunt_tb_onr(1,u)=((a(u,1))/((a(u,1))-1));
        b(tb(u),fb(u))=(1/((shunt_tb_onr(1,u))*(z(fb(u),tb(u))))); 
        z(fb(u),tb(u))=z(fb(u),tb(u))*(a(u,1)); 
        z(tb(u),fb(u))=z(fb(u),tb(u)); 
    end
end

% Including the effect(s) of static shunt-capacitance

for u=1:nbus
    if(B_shunt(u)~=0.0)
        sum(1,u)=B_shunt(u);
    end
end
for u=1:nbus
    for j=1:nbus
      if(z(u,j)==0.0)
          Ybus(u,j)=0.0;       else
          Ybus(u,j)=-(1/z(u,j)); 
      end
    end
end
for u=1:nbus
    for j=1:nbus
      sum(1,u)=sum(1,u)+b(u,j)-Ybus(u,j);
      if(j==nbus)
          Ybus(u,u)=sum(1,u); 
      end
    end
end
Gbus=real(Ybus); 
% Generating the susceptance matrix[B] by taking the imaginary part of the [Ybus] element(s)
Bbus=imag((+1)*(Ybus)); 
for u=1:(nbus-1)
    for j=1:(nbus-1)
        BforP(u,j)= -(Bbus((u+1),(j+1)));
    end
    
end
for u=1:(nbus-1)
    for j=1:(nbus-1)
        if(((bus_type(u+1))~=2))
            if(((bus_type(j+1))~=2))
               BforQ(I_count_for_BQ,K_count_for_BQ)=-(Bbus((u+1),(j+1)));
               K_count_for_BQ=K_count_for_BQ+1;
            end
        end
    end
    K_count_for_BQ=1;
    if((bus_type(u+1)~=2))
        I_count_for_BQ=I_count_for_BQ+1;
    end
end  

 fprintf('\n The updated state variables Corresponding to the iterations \n');
fprintf('Iteration number           Theta2                    Theta3                    Volatge 3\n');
fprintf('-------------------------------------------------------------------------------------------\n');
% Resetting the 'mismatch_active_vector_element_count' 
while(tolerence_checker>0)   
  mismatch_vector_element_count=0;
  mismatch_active_vector_element_count=1;
  mismatch_reactive_vector_element_count=1; 
  P_injected_sum=zeros(nbus,1);
  Q_injected_sum=zeros(nbus,1);
  tolerence_checker=0; 
P_injected_scheduled=((P_gen-P_load)/(MVA_base));
Q_injected_scheduled=((Q_gen-Q_load)/(MVA_base)); 
for u=1:nbus
    for j=1:nbus
        if((u~=j)&&(abs(Ybus(u,j)))~=0.0)
            P_injected_sum(u,1)=((P_injected_sum(u,1))+(((abs((V_new(u,1))*(V_new(j,1))*(Ybus(u,j))))*(cos((angle(Ybus(u,j)))+(angle(V_new(j,1)))-(angle(V_new(u,1)))))))); % calculating the mutual sumation part
        end
        if(j==nbus)
            P_injected_bus(u,1)=((((abs(V_new(u,1)))^(2))*(Gbus(u,u)))+(P_injected_sum(u,1))); % calculating the injecetd 'P' at bus-(u)
        end
       
    end
end

for u=1:nbus
    if(bus_type(u)~=1)
        inj_active_pow_mismatch_vector(mismatch_active_vector_element_count,1)=((P_injected_scheduled(u,1))-(P_injected_bus(u,1)))/abs(V_new(u,1)); % calculating the mismatch vector for active power injecttion
        mismatch_active_vector_element_count=mismatch_active_vector_element_count+1; 
    end
end

for u=1:(nbus-1)
    correction_voltage_angle(u,1)=0.0; 
end
correction_voltage_angle=(inv(BforP))*(inj_active_pow_mismatch_vector); 

V_new_for_Q=V_new;

for u=1:(nbus-1)
    voltage_angle_updated=angle((V_new(u+1,1)))+(correction_voltage_angle(u,1)); 
   
    V_new_for_Q(u+1,1)=abs(V_new(u+1,1))*((cos(voltage_angle_updated))+(sin(voltage_angle_updated)*(1i)));
end

for u=1:nbus
    for j=1:nbus
       
        if((u~=j)&&((abs(Ybus(u,j)))~=0.0)&&(bus_type(u)~=2))
            Q_injected_sum(u,1)=((Q_injected_sum(u,1))+(((abs((V_new(u,1))*(V_new(j,1))*(Ybus(u,j))))*(sin((angle(Ybus(u,j)))+(angle(V_new_for_Q(j,1)))-(angle(V_new_for_Q(u,1)))))))); % calculating the mutual sumation part
        end
      
        if((j==nbus)&&(bus_type(u)~=2))
            Q_injected_bus(u,1)=(-((((abs(V_new(u,1)))^(2))*(Bbus(u,u)))+(Q_injected_sum(u,1)))); % calculating the injecetd 'Q' at bus-(u)
        end
    end
end
for u=1:nbus

    if((bus_type(u)~=1)&&(bus_type(u)~=2))
        inj_reactive_pow_mismatch_vector(mismatch_reactive_vector_element_count,1)=((Q_injected_scheduled(u,1))-(Q_injected_bus(u,1)))/abs(V_new(u,1)); % calculating the mismatch vector for reactive power injecttion
        mismatch_reactive_vector_element_count=mismatch_reactive_vector_element_count+1; % updating the value of 'mismatch_vector_element_count'
    end
end

for u=1:(nbus-1)
    if((abs(inj_active_pow_mismatch_vector(u,1)))>(tolerence))
        tolerence_checker=tolerence_checker+1;
    end
end
for u=1:(nbus-1-nPV)
    if((abs(inj_reactive_pow_mismatch_vector(u,1)))>(tolerence))
         tolerence_checker=tolerence_checker+1;
    end
        
end
for u=1:(nbus-1-nPV)
    correction_voltage_magnitude(u,1)=0.0;
end
correction_voltage_magnitude=(inv(BforQ))*(inj_reactive_pow_mismatch_vector); 

correction_voltage_magnitude_count=1; 
for u=1:(nbus-1)
   
    if(bus_type(u+1)~=2)
        voltage_magnitude_updated=abs((V_new(u+1,1)))+(correction_voltage_magnitude(correction_voltage_magnitude_count,1)); % calculating the updated voltage-magnitude
        correction_voltage_magnitude_count=correction_voltage_magnitude_count+1;
    end
    if(bus_type(u+1)~=2)
        V_new(u+1,1)=voltage_magnitude_updated*((cos(angle(V_new_for_Q(u+1,1))))+(sin(angle(V_new_for_Q(u+1,1)))*(1i))); % calculating the updated voltage
    end
    if(bus_type(u+1)==2)
        V_new(u+1,1)=abs(V_new(u+1,1))*((cos(angle(V_new_for_Q(u+1,1))))+(sin(angle(V_new_for_Q(u+1,1)))*(1i))); % calculating the updated voltage
    end
end

iteration=iteration+1; 

% updating the value of iteration-count by (1)

fprintf('      %d                (%f)rad.              (%f)rad.             (%f)\n' ,iteration,((angle(V_new(2,1)))),((angle(V_new(3,1)))), abs(V_new(3,1)));
end

for u=1:nbus
   for j=1:nbus
       if((u~=j)&&(Ybus(u,j))~=0.0)
           complex_flow_line(1,flow_count)=((conj(V_new(u,1))*((V_new(u,1))-(V_new(j,1))))*(-Ybus(u,j)))+(conj(V_new(u,1))*(V_new(u,1))*(b(u,j))); % Calculating the line-flows
           line_flows(u,j)=conj(complex_flow_line(1,flow_count));
           active_flow_line(1,flow_count)=real(complex_flow_line(1,flow_count)); 
           reactive_flow_line(1,flow_count)=(-(imag(complex_flow_line(1,flow_count)))); 
          
           bus_power_injection(1,u)=bus_power_injection(1,u)+conj(complex_flow_line(1,flow_count));
           bus_power_mismatch(1,u)=((S_gen(u,1)-S_load(u,1))/MVA_base)-(bus_power_injection(1,u)); 
           flow_count=flow_count+1; 
       end
   end
end


fprintf('\n The bus-power injection(s) are :\n');
fprintf('Corresponding Real Power & Reactive power \n')
fprintf('Bus-code        Real Power(p.u.)   Reactive Power(p.u.)\n');
fprintf('--------------------------------------------------------\n');
% Displaying the bus-power injections
for u=1:nbus
    fprintf('  %d             (%f)            j(%f)\n',u,real(bus_power_injection(u)),imag(bus_power_injection(u))); 
end
