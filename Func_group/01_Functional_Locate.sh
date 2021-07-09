for file in *xyz;do
    echo \#Substrate ${file%%.*} >> func_temp
    python 01_Functional_Locator.py $file >> func_temp
    mv func_temp ${file%%.*}.txt
done
