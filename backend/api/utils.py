import nltk

nltk.download("punkt_tab")
nltk.download("averaged_perceptron_tagger_eng")


def get_nouns(text: str) -> list[str]:
    tokens = nltk.word_tokenize(text)
    tagged = nltk.pos_tag(tokens)
    nouns = [word for word, pos in tagged if pos.startswith("NN")]
    return nouns
