import streamlit as st
import pandas as pd
import numpy as np
import pickle
import joblib
from tensorflow.keras.models import load_model

model = load_model('model.h5')

def tokenization(seqs, mx_len):
	from keras_preprocessing import text, sequence
	from keras.preprocessing.text import Tokenizer

	#create and fit tokenizer
	tokenizer = Tokenizer(char_level=True)
	tokenizer.fit_on_texts(seqs)
	#represent input data as word rank number sequences
	X = tokenizer.texts_to_sequences(seqs)
	X = sequence.pad_sequences(X, maxlen=mx_len)

	X = pd.DataFrame(X)

	return X

def preprocessing(df):
	mod_s = []
	w = 7
	mx_len = w*2 + 1
	c = 0
	for i in range(len(df)):
	  x = df['Sequence'][i]
	  y = int(df['Residue No.'][i])
	  r = int(y+w) if ((y+w)<=len(x)) else (len(x))
	  l = int(y-w-1) if ((y-w-1)>=0) else 0
	  # print((l, r))
	  c+=1
	  temp =  x[l:r]
	  if(l==0 and r==(len(x))):
	    temp = '-'*(w+1-y)+temp+'-'*(y+w-len(x))
	  elif(r==(len(x))):
	    temp = temp+'-'*(y+w-len(x))
	  elif(l==0):
	    temp = '-'*(w+1-y)+temp
	  if(len(temp)!=mx_len):
	    print("Check ", temp, c)
	  mod_s.append(temp)
	# len(mod_s)
	df['Mod Sequence'] = mod_s
	df = df[['EntryID',	'Sequence', 'Mod Sequence',	'Residue No.','Modification']]
	df1 = df.drop(['EntryID', 'Sequence', 'Modification'], axis=1)

	seqs = list(df1['Mod Sequence'].values)
	X = tokenization(seqs, mx_len)

	return df1,X

def main():
	st.title("Histify")
	st.write("### Histidine Multi-class Predictions Model")

	option = st.radio('Choose Input Type', ['Single Sequence Prediction', 'Multiple Sequence Prediction'])

	if option=="Single Sequence Prediction":
		seq = st.text_input('Enter protien Sequence')
		res = st.number_input('Enter Residue number for Hustidine', step=1, value=0, format="%d")
		if st.button("Process"):
			data = [[1, seq, res, 'NA']]
			df = pd.DataFrame(data, columns=['EntryID', 'Sequence', 'Residue No.', 'Modification'])
			st.write(df)
			st.write("After preprocessing...")
			df1 = preprocessing(df)[0]

			# Display the DataFrame in the Streamlit app
			st.write(df1)

			X = preprocessing(df)[1]

			y_pred = model.predict(X)

			pred = pd.DataFrame(y_pred)
			pred.columns = ['Acetylation', 'Ribosylation', 'glycoprotein', 'hydroxylation', 'methylation', 'oxidation', 'pHOSPORYLATON', 'protein-splicing']
			pred['Prediction'] = pred.idxmax(axis=1)

			final = df1
			final['Prediction'] = pred['Prediction']
			st.write(final)


	else:
		uploaded_file = st.file_uploader('Upload a CSV', type='csv')

		if uploaded_file is not None:
			data = pd.read_csv(uploaded_file)

			st.write(data)

			st.write("After preprocessing...")
			df1 = preprocessing(data)[0]

			# Display the DataFrame in the Streamlit app
			st.write(df1)

			X = preprocessing(data)[1]

			y_pred = model.predict(X)

			y_test = np.zeros(list(y_pred.shape), dtype = int)

			for i in range(len(y_test)):
			  for j in range(8):
			      if j==6:
			        y_test[i][j]=1

			score = model.evaluate(X, y_test, verbose=0)

			st.write("%s: %.2f%%" % (model.metrics_names[1], score[1]*100))

			pred = pd.DataFrame(y_pred)
			pred.columns = ['Acetylation', 'Ribosylation', 'glycoprotein', 'hydroxylation', 'methylation', 'oxidation', 'pHOSPORYLATON', 'protein-splicing']
			pred['Prediction'] = pred.idxmax(axis=1)

			final = df1
			final['Prediction'] = pred['Prediction']
			final["EntryID"] = data["EntryID"]
			final["Modification"] = data['Modification']
			final = final[["EntryID", "Mod Sequence", "Modification", "Prediction"]]

			st.write(final)


if __name__=="__main__":
    main()