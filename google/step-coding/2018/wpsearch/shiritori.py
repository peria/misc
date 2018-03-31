#!/usr/bin/python3

import wp
import json


class Shiritori(object):

    def __init__(self):
        self.collection = wp.WikipediaCollection("data/wp.db")
        self.index = wp.Index("data/index.db", self.collection)
        self.category = None
        self.initial = None
        self.used = set()

    def reset(self):
        self.category = None
        self.initial = None
        self.used = set()

    def initialize(self, query):
        if self.index.search_category(query):
            self.category = query
            print('Category:{0}'.format(self.category))
            return self.reply('{0} のカテゴリでしりとりしましょう。何から初めますか。'.format(query))
        else:
            return self.reply('そのカテゴリは知りません。'.format(query))

    def answer(self, query):
        document = self.collection.get_document_by_id(query)
        if not document:
            return False, self.reply('そんな単語はありません。私の勝ちです')
        word = self.index.get_word(self.category, query)
        if not word:
            return False, self.reply('{0}は{1}ではありません。私の勝ちです'.format(
                query, self.category))
        if self.initial and word.reading[0] != self.initial:
            return False, self.reply('{0}は {1} から始まりません。私の勝ちです'.format(
                query, self.initial))
        if query in self.used:
            return False, self.reply('既に {0}は言われています。私の勝ちです'.format(
                query))
        self.used.add(query)

        # Accept user's answer
        last = query[-1]
        words = self.index.search(self.category, last)
        for word in words:
            if word.id in self.used:
                continue
            self.used.add(word.id)
            return True, self.reply(word.id)

        self.reset()
        return False, self.reply('もうありません。私の負けです')

    def reply(self, query='よく分かりません'):
        return json.dumps({
            'textToSpeech': query
        }, indent=2, separators=(',', ': '), ensure_ascii=False)

    @property
    def is_initialized(self):
        return bool(self.category)
