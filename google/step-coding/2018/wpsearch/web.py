import bottle
import json
import os
from shiritori import Shiritori

# Statefull
g_shiritori = Shiritori()

@bottle.route('/action')
def action():
   global g_shiritori
   bottle.response.content_type = 'application/json'

   query = bottle.request.query.q
   # Easter eggs
   # Thanks, goodness
   if query == 'ありがとう':
      return g_shiritori.thanks()
   if query == 'カテゴリー' or query == 'カテゴリ':
      return g_shiritori.ask_categories()

   # Set category
   if not g_shiritori.is_initialized:
      return g_shiritori.initialize(query)

   is_cont, reply = g_shiritori.answer(query)
   if not is_cont:
      g_shiritori.reset()
   return reply


@bottle.route('/article/<title>')
def article(title):
    article = collection.get_document_by_id(title)
    bottle.response.content_type = 'application/json'
    if article is None:
        bottle.abort(404, "Not found")
    return json.dumps({
        'title': article.title,
        'text': "<<<Omitted>>>",
        'opening_text': article.opening_text,
        'auxiliary_text': article.auxiliary_text,
        'categories': article.categories,
        'headings': article.headings,
        'wiki_text': "<<<Omitted>>>",
        'popularity_score': article.popularity_score,
        'num_incoming_links': article.num_incoming_links,
    }, indent=2, separators=(',', ': '), ensure_ascii=False)

@bottle.route('/article/wiki_text/<title>')
def article_wiki_text(title):
    article = collection.get_document_by_id(title)
    if article is None:
        bottle.abort(404, "Not found")
    bottle.response.content_type = 'text/plain;charset=utf-8'
    return article.wiki_text

@bottle.route('/article/text/<title>')
def article_text(title):
    article = collection.get_document_by_id(title)
    if article is None:
        bottle.abort(404, "Not found")
    bottle.response.content_type = 'text/plain;charset=utf-8'
    return article.text()

port = 8081
if 'WPSEARCH_PORT' in os.environ:
    port = int(os.environ['WPSEARCH_PORT'])
bottle.run(host='0.0.0.0', port=port, reloader=False)
