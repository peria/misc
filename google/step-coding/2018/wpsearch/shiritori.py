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
        self.limit = 10

    def reset(self):
        self.category = None
        self.initial = None
        self.used = set()
        self.limit = 10

    def initialize(self, query):
        words = self.index.search_category(query)
        if words:
            self.category = query
            print('Category:{0}'.format(self.category))
            # Output hints
            print([word.id for word in words])
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
        words = self.index.search(self.category, last, self.limit)
        for word in words:
            if word.id in self.used:
                continue
            self.used.add(word.id)
            self.initial = word.reading[-1]
            return True, self.reply(word.id)

        return False, self.reply('もうありません。私の負けです')

    def thanks(self):
        if self.limit > 2:
            self.limit -= 1
        print('LIMIT {0}'.format(self.limit))
        return self.reply('どういたしまして')

    def reply(self, query='よく分かりません'):
        return json.dumps({
            'textToSpeech': query
        }, indent=2, separators=(',', ': '), ensure_ascii=False)

    @property
    def is_initialized(self):
        return bool(self.category)
