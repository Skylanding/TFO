## TFO Model Document

![TFO model](D:\Desktop\TFO_NOTES\Document\TFO model.png)



The whole model framework can be divided into 4 sections:

1. Construct 3D vessel model, a dimension is used to simulate the time axis.

2. Obtain TOF curve with Monte Carlo simulation.

3. Use the light intensity calculated by TOF to draw the beating signal curve.

4. Substitute the peak-valley value into BLL to calculate SpO2.









#### SpO2 in one subject

$$
SpO2 = \frac{[HbO2]}{[HbO2]+[HbO]} = \frac{c(HbO2)}{c(HbO2)+c(HbO)} = \frac{ratio}{ratio+1}(1)
$$

$$
S p O 2=\frac{\epsilon_{H b, \lambda_{1}}-R * \epsilon_{H b, \lambda_{2}}}{\epsilon_{H b, \lambda_{1}}-\epsilon_{H b O_{2}, \lambda_{1}}+R *\left(\epsilon_{H b O_{2}, \lambda_{2}}-\epsilon_{H b, \lambda_{2}}\right)}
(2)
$$



It can be seen from formula (2) that Spo2 is related to R value. And about R value:
$$
R = \frac{ΔAλ1}{ΔAλ2}\
 = \frac{\sum_jΔμa_j(λ1)*<L_j>}{\sum_jΔμa_j(λ2)*<L_j>} (Expansion 1,\ j\ means\ different \ layers)\\
 = \frac{\sum_i(Δc_i*\epsilon_{i}*L)_{λ1}}{\sum_i(Δc_i*\epsilon_{i}*L)_{λ2}} (Expansion 2,\ i\ means\ different\ chromophores)
$$



##### For expansion 1:

$$
R = \frac{Δμa_1(λ1)*<L_1>_{λ1}}{Δμa_1(λ2)*<L_1>_{λ2}}(Consider\ only\ Maternal\ Layer\ changes)\\
  = \frac{μa_1(λ1)*<L_1>_{λ1}}{μa_1(λ2)*<L_1>_{λ2}}(The\ change\ is\ the\ same\ at\ both\ wavelengths)\\
  = \frac{(L_1)_{λ1}}{(L_1)_{λ2}} * \frac{μa_1(λ1)}{μa_1(λ2)}
$$
That's the reverse process in forward model. Because in general situations we set up different absorption coefficient across different layers. This formula supports our idea in step1 and step2 of one-subject model.



##### For expansion 2:  

(c1: concentration of Hb c2: concentration of HbO2)
$$
R = \frac{\sum_i(Δc_i*\epsilon_{i}*L)_{λ1}}{\sum_i(Δc_i*\epsilon_{i}*L)_{λ2}}\\
  = \frac{\sum(Δc_1*\epsilon_{1}*L)_{λ1}+\sum(Δc_2*\epsilon_{2}*L)_{λ1}}{\sum(Δc_1*\epsilon_{1}*L)_{λ2}+\sum(Δc_2*\epsilon_{2}*L)_{λ2}}\\
  (Consider\ Layer1\ change)\\
  = \frac{\sum(c_1*\epsilon_{1}*L)_{λ1}+\sum(c_2*\epsilon_{2}*L)_{λ1}}{\sum(c_1*\epsilon_{1}*L)_{λ2}+\sum(c_2*\epsilon_{2}*L)_{λ2}}\\
  = \frac{\sum(L_{λ1}*(\epsilon_{1}+ratio*\epsilon_{2})_{λ1})}{\sum(L_{λ2}*(\epsilon_{1}+ratio*\epsilon_{2})_{λ2})}\\
  = \frac{\sum(\epsilon_{1,λ1}+ratio*\epsilon_{2,λ1})}{\sum(\epsilon_{1,λ2}+ratio*\epsilon_{2,λ2})}*(\frac{L_{λ1}}{L_{λ2}})
$$



The Lppath values of the two expansions are measured based on the Monte Carlo simulation of the following given data.

```
# (1) 735nm
# Maternal Abdominal Wall: mu_a=0.0094(mm^-1) mu_s=13.22(mm^-1)
# Maternal Uterus: mu_a=0.016 (mm^-1) mu_s=10.8(mm^-1)
# Amniotic Fluid: mu_a=0.0025 (mm^-1) mu_s=0.1(mm^-1)
# Fetal Tissues: mu_a=0.0187 (mm^-1) mu_s=12.33(mm^-1)

# (1) 850nm
# Maternal Abdominal Wall: mu_a=0.009(mm^-1) mu_s=12(mm^-1)
# Maternal Uterus: mu_a=0.01 (mm^-1) mu_s=8.15(mm^-1)
# Amniotic Fluid: mu_a=0.0042 (mm^-1) mu_s=0.1(mm^-1)
# Fetal Tissues: mu_a=0.013 (mm^-1) mu_s=9.916(mm^-1)
```

After checking and validation, only less than 0.7% error between ratio in formula (1) and formula (2), which supports that Spo2 depends on ratio. So if we wan t to check the connection between Spo2 and μa, we can focus on **ratio**.





##### About μa, for overall μa:

$$
μa = \epsilon_{Hb}*[Hb]+\epsilon_{HbO2}*[HbO2]\\
   = β*[\epsilon_{Hb}*c(Hb)+\epsilon_{HbO2}*c(HbO2)]\ (β\ is\ volume\ factor)\\
   = β*(\epsilon_{Hb}+\epsilon_{HbO2}*ratio)*c(Hb)
$$
In this formula, the value of μa depends on β, ratio and the concentration of Hb. We have proved Spo2 depends on ratio, so we need to get β and c(Hb) to establish clear relationship between μa and Spo2.

In section 1 we have constructed a vessel model to simulate changes in the state of blood vessel in one cycle. By scanning the cross-sectional area of the vessel from the starting position, the curve can be obtained. 

<img src="D:\Desktop\TFO_NOTES\fig\output.png" alt="output" style="zoom:67%;" />



About c(Hb), the concentration of Hb for specific group is in a range, through searching around and make verification with our collaborator's data, we can determine a range for analysis, this range is 12-16.



<img src="C:\Users\22877\AppData\Roaming\Typora\typora-user-images\image-20230906103719181.png" alt="image-20230906103719181" style="zoom:50%;" />



After determining the range, we can perform validation on several values within the range, each concentration value representing an individual. The value of spo2 is fixed at 0.98, and the change curve of μa can be fitted by substituting the formula into the calculation.

<img src="C:\Users\22877\AppData\Roaming\Typora\typora-user-images\image-20230906104356100.png" alt="image-20230906104356100" style="zoom: 25%;" />
