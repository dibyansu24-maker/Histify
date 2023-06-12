import streamlit as st
import pandas as pd
import numpy as np
import pickle
import joblib
import streamlit.components.v1 as components
from tensorflow.keras.models import load_model

st.set_page_config(
	page_title="Histify", 
	page_icon="🧬", 
	layout="wide",
	menu_items={
        'Get Help': 'https://github.com/dibyansu24-maker/Histify',
        'Report a bug': "https://github.com/dibyansu24-maker/Histify/issues",
        'About': "## Sequence-based multiple histidine function prediction methods trained on deep neural network"
    }
    )

@st.cache_resource
def loaded_model():
	return load_model('model.h5')

model = loaded_model()

@st.cache_data
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

@st.cache_data
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
	df = df[['EntryID',	'Sequence', 'Mod Sequence',	'Residue No.']]
	df1 = df.drop(['EntryID', 'Sequence'], axis=1)

	seqs = list(df1['Mod Sequence'].values)
	X = tokenization(seqs, mx_len)

	return df1,X

@st.cache_data
def SingleSequencePred(seq, res):
	data = [[1, seq, res]]
	df = pd.DataFrame(data, columns=['EntryID', 'Sequence', 'Residue No.'])
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

	st.write('**Predictions:**')
	final = df1
	final['Prediction'] = pred['Prediction']
	st.write(final)

@st.cache_data
def MultiSequencePred(data):
	st.write("After preprocessing...")
	df1 = preprocessing(data)[0]

	# Display the DataFrame in the Streamlit app
	st.write(df1)

	X = preprocessing(data)[1]

	y_pred = model.predict(X)

	pred = pd.DataFrame(y_pred)
	pred.columns = ['Acetylation', 'Ribosylation', 'glycoprotein', 'hydroxylation', 'methylation', 'oxidation', 'pHOSPORYLATON', 'protein-splicing']
	pred['Prediction'] = pred.idxmax(axis=1)

	st.write('**Predictions:**')
	final = df1
	final['Prediction'] = pred['Prediction']
	final["EntryID"] = data["EntryID"]
	final = final[["EntryID", "Mod Sequence", "Prediction"]]

	st.write(final)

def main():
	st.title("Histify")
	st.write("### Histidine Multi-class Predictions Model")
	st.caption('Sequence-based multiple histidine function prediction methods trained on deep neural network')
	st.divider()

	tab1, tab2 = st.tabs(['Single Sequence Prediction', 'Multiple Sequence Prediction'])

	# option = st.radio('Choose Input Type', ['Single Sequence Prediction', 'Multiple Sequence Prediction'])

	with tab1:
		with st.expander("See Example", expanded=True):
			st.write("*Example Sequence: MV**H**LTPEEKSAVTAL*")
			st.write("*Example Residue No.: 3*")

			st.markdown('<p style="color:grey;">Definition:</p>', unsafe_allow_html=True)
			st.markdown("<ul style='color:grey;'><li>'Sequence' represents protein sequence.</li> <li>'Residue No.' represents the position where histidine is present.</li></ul>", unsafe_allow_html=True)		    
		
		seq = st.text_input('**Enter protein Sequence**')
		res = st.number_input('**Enter Residue number for Histidine**', step=1, value=0, format="%d")
		st.info("To ensure optimal performance, we have employed a window size of 7 in our model.", icon="ℹ️")

		
		if st.button("Process"):
			SingleSequencePred(seq, res)
			st.success('Prediction done successfully!', icon="✅")
		if st.button("Reset"):
			st.experimental_rerun()


	with tab2:
		# Example dataframe in CSV format
		example_dataframe = pd.DataFrame({
		    'Sequence': ['MSIPETQKGVIFYESHGKLEYKDI', 'MVHLTPEEKSAVTAL', 'MVLSPADKTNVKAAWGKVGAHAGEY'],
		    'Residue No.': [16, 3, 21]
		})
		example_dataframe = example_dataframe.style.set_properties(**{'white-space': 'nowrap'})
		
		with st.expander("See Example", expanded=True):
			st.write("*Example Dataframe:*")
			st.write(example_dataframe)
			st.markdown('<p style="color:grey;">Definition:</p>', unsafe_allow_html=True)
			st.markdown("<ul style='color:grey;'><li>'Sequence' represents protein sequence.</li> <li>'Residue No.' represents the position where histidine is present.</li></ul>", unsafe_allow_html=True)

		uploaded_file = st.file_uploader('**Upload a CSV**', type='csv')
		st.info("To ensure optimal performance, we have employed a window size of 7 in our model.", icon="ℹ️")

		# Display the example dataframe if no file uploaded yet
		if uploaded_file is None:
			st.warning("Please use 'Sequence' and 'Residue No.' as column names in the CSV file to avoid errors.", icon="⚠️")
		else:
			data = pd.read_csv(uploaded_file)
			st.session_state["uploaded_file"] = data
			if 'EntryID' not in data.columns:
				data.insert(loc=0, column='EntryID', value=np.arange(data.shape[0]))

			data = data[['EntryID', 'Sequence', 'Residue No.']]

			data_styl = data.style.set_properties(**{'white-space': 'nowrap'})
			st.write(data_styl)

			MultiSequencePred(data)
			st.success('Prediction done successfully!', icon="✅")

	st.markdown("[Check our complete Github repo here.](https://github.com/dibyansu24-maker/Histify)")


if __name__=="__main__":
    main()